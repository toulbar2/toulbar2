import csv
import sys


page_title = "Results"
default_cutoff = 3600

# Header constants from the results file
bestsol_postfix = "_bestsol"
cputime_postfix = "_cputime"
status_postfix = "_status"
filename_column = "Problem"
probstat_columns = ["nbvar", "max_dom", "nbconstr", "max_arity"]

# CSS classes
timed_out_class = "timed_out"
opt_best_class = "opt_best"
opt_ok_class = "opt_ok"
err_class = "error"


class SolverResult(object):

    def __init__(self, bestsol, cputime, status):
        self.bestsol_text = bestsol
        if bestsol == "?":
            self.bestsol = sys.maxint
        else:
            self.bestsol = int(bestsol)

        self.cputime_text = cputime  # Keep a copy of the actual text for display later
        try:
            self.cputime = float(cputime)
            self.cputime_text = str(int(round(self.cputime)))
        except:
            self.cputime = default_cutoff

        self.status = status

    def timed_out(self):
        return self.status == "UNK"

    def is_opt(self):
        return self.status == "SC"

    def is_err(self):
        return self.cputime_text in ["MZN", "mem", "ARITY", "32-bit"]

    def get_ratio_to_best(self, best):
        try:
            if self.bestsol != 0:
                return (1.0 - (float(self.bestsol) - float(best.bestsol)) / float(self.bestsol)) * 0.5
            else:
                return 1.0

        except Exception as e:
            print self.status, best.status, self.bestsol, best.bestsol
            raise e

    def __lt__(self, other):
        # Used for sorting results
        if (other.status == "UNK" and self.status[0] == "S") or\
           (other.status == "S" and self.status == "SC"):
            return True

        if self.bestsol < other.bestsol:
            return True
        
        if self.bestsol == other.bestsol and self.cputime < other.cputime:
            return True

        return False

    def __str__(self):
        return "%d,%s,%s" % (self.bestsol, self.cputime_text, self.status)


def read_results(filename):
    solver_names = []
    rows = []  # Each row is a tuple of four items (problemfilename, dictionary of problem stats, dictionary of solver name to SolverResult objects, best SolverResult)
    with open(filename, "rt") as f:
        reader = csv.DictReader(f)

        # Get the list of solver names
        fieldnames = reader.fieldnames
        for x in fieldnames:
            if x != filename_column and x not in probstat_columns:
                solver_name = x
                solver_name = solver_name.replace(bestsol_postfix, "")
                solver_name = solver_name.replace(cputime_postfix, "")
                solver_name = solver_name.replace(status_postfix, "")
                if solver_name not in solver_names:
                    solver_names.append(solver_name)

        for row in reader:
            problem = row[filename_column]
            stats = dict((name, int(row[name])) for name in probstat_columns)

            row_results = {}
            for solver_name in solver_names:
                bestsol = row["%s%s" % (solver_name, bestsol_postfix)]
                cputime = row["%s%s" % (solver_name, cputime_postfix)]
                status = row["%s%s" % (solver_name, status_postfix)]

                res = SolverResult(bestsol, cputime, status)
                row_results[solver_name] = res

            best_result = min(row_results.values())
            rows.append((problem, stats, row_results, best_result))

    return solver_names, rows


# Returns a HTML string of the problem stats
def format_stats(stats):
    ret = "<table class='stat'>"
    ret += "".join("<tr><td class='title'>%s:</td> <td class='value'>%d</td></tr>" % (p, stats[p]) for p in probstat_columns)
    ret += "</table>"
    return ret


# Returns a HTML string for a single result cell
def format_result(res):
    return """<table class="result">
<tr><td class='title'>S:</td><td class='value'>%s</td></tr>
<tr><td class='title'>O:</td><td class='value'>%s</td></tr>
<tr><td class='title'>T:</td><td class='value'>%s</td></tr>
</table>""" % (res.status, res.bestsol_text, res.cputime_text)


def html_legend(shading):
    shading_html = ""
    if shading:
        shading_html = """<tr><td style="background-color: rgba(176,255,0,0.5)">Text</td><td>Non-optimal solutions are shaded by ratio to the best.</td></tr>"""

    return """
<table class="resultstable">
<tr><th class='title'>Key</th><th>Meaning</th></tr>
<tr><td class='title'>S</td><td>
SC: Optimal solution found and proved<br/>
S: Satisfiable, solution found<br/>
UNK: No solution returned
</td></tr>
<tr><td class='title'>O</td><td>Best objective value found</td></tr>
<tr><td class='title'>T</td><td>
Time in seconds<br/>
mem: over the memory limit<br/>
MZN: error during MiniZinc flattening<br/>
ARITY: maximum arity too large<br/>
32-bit: maximum sum of costs too large<br/>
N/A: result non available
</td></tr>
</table><br>

<table class="resultstable">
<tr><td class="%s">Text</td><td>Optimal solution with the best CPU time</td></tr>
<tr><td class="%s">Text</td><td>Optimal solution within time limit</td></tr>
%s
<tr><td class="%s">Text</td><td>Time out</td></tr>
<tr><td class="%s">Text</td><td>Error</td></tr>
</table>
<br/>
""" % (opt_best_class, opt_ok_class, shading_html, timed_out_class, err_class)


def html_header():
    return """<html>
<head>
<title>%s</title>
<style type="text/css">
.resultstable, .resultstable th, .resultstable td {
    border: 1px solid black;
}

.resultstable td {
    padding: 0 0.2em;
}

.stat, .result, .stat td, .result td {
    border: 0;
    border-spacing: 0;
    padding: 0;
    width: 100%%;
}

.stat .value, .result .value {
    text-align: right;
}

.%s {
    background-color: #FFDF7F;
}

.%s {
    background-color: #00FF00;
}

.%s {
    background-color: #AFFF00;
}

.%s {
    background-color: #FFAFAF;
}
</style>
</head>
<body>""" % (page_title, timed_out_class, opt_best_class, opt_ok_class, err_class)

def html_footer():
    return "</body></html>"


def format_html(solver_names, rows, out_f=sys.stdout, header_interval=sys.maxint, shading=True):
    # Outputs a formated HTML page to 'out_f'. Reprints the header row each 'header_interval' number of rows.
    # 'shading' turns on/off colouring the background in shades of ratio to the best solution.

    solver_headers = "".join("<th>%s</th>" % s for s in solver_names)
    # solver_headers = ""
    table_header = "<table class=\"resultstable\">"
    table_header_row = "<tr><th>%s</th> <th>Stats</th> <th>Best Solver</th> %s </tr>" % (filename_column, solver_headers)
    table_footer = "</table>"

    print >> out_f, html_header()
    print >> out_f, html_legend(shading)
    print >> out_f, table_header

    for i, row in enumerate(rows):
        problem, stats, row_results, best_result = row
        if i % header_interval == 0:
            print >> out_f, table_header_row
        print >> out_f, "<tr><td>%s</td><td>%s</td><td>%s</td>" % (problem, format_stats(stats), format_result(best_result))

        for solver_name in solver_names:
            res = row_results[solver_name]
            tdclass = ""
            if res.is_err():
                tdclass = err_class
            elif res.timed_out():
                tdclass = timed_out_class
            elif res is best_result:
                tdclass = opt_best_class
            elif res.is_opt():
                tdclass = opt_ok_class

            if shading and not tdclass and res.status[0] == "S":
                d = res.get_ratio_to_best(best_result)
                tdclass = "\" style=\"background-color: rgba(176,255,0,%.4f);" % d

            print >> out_f, "<td class=\"%s\">%s</td>" % (tdclass, format_result(res))

        print >> out_f, "</tr>"

    print >> out_f, table_footer
    print >> out_f, html_footer()


def main():
    solver_names, rows = read_results("results.csv")
    format_html(solver_names, rows, out_f=open("results.html", "wt"), header_interval=50)
    format_html(solver_names, rows, out_f=open("results_formatted_noshading.html", "wt"), header_interval=50, shading=False)


if __name__ == '__main__':
    main()
