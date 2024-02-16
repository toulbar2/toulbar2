/* Copyright (C) 2018-2023 INRAE @ Gauthier Quesnel
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include "lp-parser.hpp"

#include <algorithm>
#include <cassert>
#include <climits>
#include <cstring>
#include <fstream>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>

#include "core/tb2types.hpp"

extern const char* PrintFormatProb;

int baryonyx::precision;

Double baryonyx::Floor(Double v)
{
    return floorl(v + ToulBar2::epsilon);
}

Double baryonyx::Ceil(Double v)
{
    return ceill(v - ToulBar2::epsilon);
}

std::ostream& operator<<(std::ostream& os, const baryonyx::file_format_error_tag& obj)
{
    os << static_cast<std::underlying_type<baryonyx::file_format_error_tag>::type>(obj);
    return os;
}

constexpr bool are_equal(const std::string_view lhs,
    const std::string_view rhs) noexcept
{
    if (rhs.size() != lhs.size())
        return false;

    std::string_view::size_type i = 0;
    std::string_view::size_type e = lhs.size();

    for (; i != e; ++i)
        if (std::tolower(lhs[i]) != std::tolower(rhs[i]))
            return false;

    return true;
}

constexpr bool is_char_operator(const int c) noexcept
{
    return c == '<' || c == '>' || c == '=';
}

constexpr bool starts_with_operator(const std::string_view str) noexcept
{
    return (!str.empty() && is_char_operator(str[0]));
}

constexpr bool is_char_name(const int c) noexcept
{
    if ((c >= static_cast<int>('a') && c <= static_cast<int>('z')) || (c >= static_cast<int>('A') && c <= static_cast<int>('Z')) || (c >= static_cast<int>('0') && c <= static_cast<int>('9')))
        return true;

    switch (c) {
    case '!':
    case '"':
    case '#':
    case '$':
    case '%':
    case '&':
    case '(':
    case ')':
    case ',':
    case '.':
    case ';':
    case '?':
    case '@':
    case '_':
    case '{':
    case '}':
    case '~':
        return true;
    default:
        return false;
    }
}

constexpr bool starts_with_name(const std::string_view str) noexcept
{
    return (!str.empty() && is_char_name(str[0]));
}

constexpr bool is_number(const int c) noexcept
{
    return (c >= static_cast<int>('0') && c <= static_cast<int>('9')) || c == '.' || c == 'e' || c == 'E' || c == '-' || c == '+';
}

constexpr bool starts_with_a_number(const int c) noexcept
{
    return (c >= static_cast<int>('0') && c <= static_cast<int>('9')) || c == '.' || c == '-' || c == '+';
}

constexpr bool starts_with_number(const std::string_view str) noexcept
{
    if (str.empty())
        return false;

    if (str[0] == 'i' || str[0] == 'I')
        if (are_equal(str, "infinity") || are_equal(str, "inf"))
            return true;

    return is_number(str[0]);
}

constexpr bool is_separator(const int c) noexcept
{
    switch (c) {
    case '<':
    case '=':
    case '>':
    case ':':
    case '-':
    case '+':
    case '[':
    case ']':
    case '*':
    case '^':
        return true;
    default:
        return false;
    }
}

static const std::string_view keywords[] = {
    "binary", "binaries", "bin", "bound", "bounds", "general", "generals",
    "gen", "end", "st", "subject", "sush", "s.t.", "st."
};

constexpr bool is_keyword(const std::string_view str) noexcept
{
    // std::find_if is not constexpr in C++17
    // return std::find_if(std::begin(array),
    //                     std::end(array),
    //                     std::bind(are_equal, std::placeholders::_1, str)) !=
    //        std::end(array);

    for (auto i : keywords)
        if (are_equal(i, str))
            return true;

    return false;
}

struct operator_token {
    constexpr operator_token() noexcept = default;

    constexpr operator_token(baryonyx::operator_type op_, int read_) noexcept
        : op(op_)
        , read(read_)
    {
    }

    baryonyx::operator_type op = baryonyx::operator_type::equal;
    int read = 0;
};

struct real_token {
    constexpr real_token() noexcept = default;

    constexpr real_token(Double value_, int read_)
        : value(value_)
        , read(read_)
    {
    }

    Double value = 0.0;
    int read = 0;
};

struct function_element_token {
    constexpr function_element_token() noexcept = default;

    constexpr function_element_token(Double factor_, std::string_view name_,
        int read_) noexcept
        : factor(factor_)
        , name(name_)
        , read(read_)
    {
    }

    Double factor = 0.0;
    std::string_view name = {};
    int read = 0;
};

struct sub_bound_token {
    constexpr sub_bound_token() noexcept = default;

    constexpr sub_bound_token(Double value_, int read_,
        baryonyx::operator_type op_) noexcept
        : value(value_)
        , read(read_)
        , op(op_)
    {
    }

    Double value = 0.0;
    int read = 0;
    baryonyx::operator_type op = baryonyx::operator_type::equal;
};

struct bound_token {
    constexpr bound_token() noexcept = default;

    constexpr bound_token(Double min_, Double max_, std::string_view name_,
        int read_) noexcept
        : min(min_)
        , max(max_)
        , name(name_)
        , read(read_)
    {
    }

    Double min = 0.0;
    Double max = std::numeric_limits<Double>::infinity();
    std::string_view name;
    int read = 0;
};

struct raw_problem_status {
    baryonyx::file_format_error_tag tag;
    unsigned int line = 0u;
    unsigned int column = 0u;

    constexpr raw_problem_status(baryonyx::file_format_error_tag tag_) noexcept
        : tag(tag_)
    {
    }

    constexpr operator bool() const noexcept
    {
        return tag == baryonyx::file_format_error_tag::success;
    }
};

class stream_buffer {
public:
    constexpr static int stream_buffer_size = 10;

    using string_view_array = std::array<std::string_view, stream_buffer_size>;

    stream_buffer(std::istream& is_)
        : is(is_)
    {
        for (int i = 0; i != stream_buffer_size; ++i)
            buffer_ptr[i] = next_token();
    }

    void pop_front()
    {
        auto tmp = buffer_ptr[0];

        buffer_ptr[0] = buffer_ptr[1];
        buffer_ptr[1] = buffer_ptr[2];
        buffer_ptr[2] = buffer_ptr[3];
        buffer_ptr[3] = buffer_ptr[4];
        buffer_ptr[4] = buffer_ptr[5];
        buffer_ptr[5] = buffer_ptr[6];
        buffer_ptr[6] = buffer_ptr[7];
        buffer_ptr[7] = buffer_ptr[8];
        buffer_ptr[8] = buffer_ptr[9];
        buffer_ptr[9] = tmp;

        do_push_back(1);
    }

    void pop_front(int size)
    {
        assert(size >= 1 && size < stream_buffer_size);

        auto tmp = buffer_ptr[0];

        buffer_ptr[0] = buffer_ptr[(size + 0) % stream_buffer_size];
        buffer_ptr[1] = buffer_ptr[(size + 1) % stream_buffer_size];
        buffer_ptr[2] = buffer_ptr[(size + 2) % stream_buffer_size];
        buffer_ptr[3] = buffer_ptr[(size + 3) % stream_buffer_size];
        buffer_ptr[4] = buffer_ptr[(size + 4) % stream_buffer_size];
        buffer_ptr[5] = buffer_ptr[(size + 5) % stream_buffer_size];
        buffer_ptr[6] = buffer_ptr[(size + 6) % stream_buffer_size];
        buffer_ptr[7] = buffer_ptr[(size + 7) % stream_buffer_size];
        buffer_ptr[8] = buffer_ptr[(size + 8) % stream_buffer_size];
        buffer_ptr[9] = tmp;

        do_push_back(size);
    }

    void print(const std::string_view msg) const
    {
        std::fprintf(stdout, "%.*s", static_cast<int>(msg.size()), msg.data());
        for (int i = 0; i < stream_buffer_size; ++i)
            std::fprintf(stdout, "[%d: (%.*s)]", i,
                static_cast<int>(buffer_ptr[i].size()),
                buffer_ptr[i].data());
        std::fprintf(stdout, "\n");
    }

    const std::array<std::string_view, stream_buffer_size>
    array() const noexcept
    {
        return buffer_ptr;
    }

    std::string_view first() const noexcept { return buffer_ptr[0]; }

    std::string_view second() const noexcept { return buffer_ptr[1]; }

    std::string_view third() const noexcept { return buffer_ptr[2]; }

    std::string_view fourth() const noexcept { return buffer_ptr[3]; }

    unsigned line() const noexcept { return 0u; }

    unsigned column() const noexcept { return 0u; }

private:
    struct stream_token {
        char buffer[512] = { '\0' };
        std::size_t current = 0;
    };

    std::array<std::string_view, stream_buffer_size> buffer_ptr;
    std::array<stream_token, stream_buffer_size> token;
    std::istream& is;
    int current_token_buffer = 0;

    void do_push_back(int size)
    {
        assert(size > 0 && size <= stream_buffer_size);

        for (int i = stream_buffer_size - size; i < stream_buffer_size; ++i)
            buffer_ptr[i] = next_token();
    }

    std::string_view next_token() noexcept
    {
        if (token[current_token_buffer]
                .buffer[token[current_token_buffer].current]
            == '\0') {
            current_token_buffer = (current_token_buffer + 1) % stream_buffer_size;
            if (!is.good())
                return std::string_view();

            token[current_token_buffer].buffer[0] = '\0';
            while (is >> token[current_token_buffer].buffer) {
                if (token[current_token_buffer].buffer[0] == '\\') {
                    is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    continue;
                }

                if (token[current_token_buffer].buffer[0] == '\0')
                    continue;

                break;
            }

            token[current_token_buffer].current = 0;
            if (token[current_token_buffer].buffer[0] == '\0')
                return std::string_view();
        }

        if (token[current_token_buffer]
                .buffer[token[current_token_buffer].current]
            == '\0')
            return std::string_view();

        bool starts_with_number = starts_with_a_number(token[current_token_buffer]
                                                           .buffer[token[current_token_buffer].current]);

        std::size_t start = token[current_token_buffer].current++;

        if (is_separator(token[current_token_buffer].buffer[start]))
            return std::string_view(&token[current_token_buffer].buffer[start], 1);

        if (token[current_token_buffer]
                .buffer[token[current_token_buffer].current]
            == '\0')
            return std::string_view(&token[current_token_buffer].buffer[start], 1);

        while (token[current_token_buffer]
                   .buffer[token[current_token_buffer].current]
            != '\0') {
            if ((!starts_with_number ||
                 (token[current_token_buffer].buffer[token[current_token_buffer].current - 1] != 'e'
                  &&
                  token[current_token_buffer].buffer[token[current_token_buffer].current - 1] != 'E'
                 )
                )
                &&
                is_separator(token[current_token_buffer].buffer[token[current_token_buffer].current]))
                break;

            if (starts_with_number && !is_number(token[current_token_buffer].buffer[token[current_token_buffer].current]))
                break;

            ++token[current_token_buffer].current;
        }

        //cout << std::string_view(&token[current_token_buffer].buffer[start],
        //    token[current_token_buffer].current - start) << endl;

        return std::string_view(&token[current_token_buffer].buffer[start],
            token[current_token_buffer].current - start);
    }
};

struct problem_parser {
    baryonyx::raw_problem m_problem;
    std::unordered_map<std::string_view, int> m_cache_variable;

    problem_parser()
    {
        m_problem.strings = std::make_shared<baryonyx::string_buffer>();
    }

    int get_or_assign_variable(const std::string_view value)
    {
        if (auto it = m_cache_variable.find(value); it != m_cache_variable.end())
            return it->second;

        if (std::cmp_greater(m_cache_variable.size(), INT_MAX))
            return -1;

        auto size = static_cast<int>(m_cache_variable.size());
        auto string = m_problem.strings->append(value);

        m_cache_variable.try_emplace(string, size);

        m_problem.vars.names.emplace_back(string);
        m_problem.vars.values.emplace_back(0, std::numeric_limits<int>::max(),
            baryonyx::variable_type::real);

        return size;
    }

    int get_variable(const std::string_view value) const noexcept
    {
        if (auto it = m_cache_variable.find(value); it != m_cache_variable.end())
            return it->second;

        return -1;
    }

    [[nodiscard]] bool append_to_objective(const function_element_token& elem)
    {
        if (elem.name.empty()) {
            m_problem.objective.value += elem.factor;
            return true;
        } else {
            auto id = get_or_assign_variable(elem.name);
            if (id == -1)
                return false;

            auto it = std::find_if(
                std::begin(m_problem.objective.elements),
                std::end(m_problem.objective.elements),
                [id](const auto& obj_elem) { return obj_elem.variable_index == id; });

            if (it == std::end(m_problem.objective.elements))
                m_problem.objective.elements.emplace_back(elem.factor, id);
            else
                it->factor += elem.factor;

            return true;
        }
    }

    [[nodiscard]] bool
    set_boolean_variable(const std::string_view value) noexcept
    {
        auto id = get_variable(value);
        if (id < 0)
            return false;

        m_problem.vars.values[id].type = baryonyx::variable_type::binary;
        m_problem.vars.values[id].min = 0;
        m_problem.vars.values[id].max = 1;

        return true;
    }

    [[nodiscard]] bool
    set_integer_variable(const std::string_view value) noexcept
    {
        auto id = get_variable(value);
        if (id < 0)
            return false;

        m_problem.vars.values[id].type = baryonyx::variable_type::general;

        return true;
    }

    [[nodiscard]] bool set_bound_variable(bound_token bound)
    {
        auto id = get_variable(bound.name);
        if (id < 0)
            return false;

        // We only use a binary/integer solver.
        // m_problem.vars.values[id].min = bound.min;
        // m_problem.vars.values[id].max = bound.max;

        m_problem.vars.values[id].min = (bound.min == -std::numeric_limits<Double>::infinity()
                ? std::numeric_limits<int>::min()
                : static_cast<int>(baryonyx::Ceil(bound.min)));

        m_problem.vars.values[id].max = (bound.max == std::numeric_limits<Double>::infinity()
                ? std::numeric_limits<int>::max()
                : static_cast<int>(baryonyx::Floor(bound.max)));

        return true;
    }
};

// constexpr copy_n problem
std::optional<Double> read_real(const std::string_view buf) noexcept
{
    if (buf.size() >= 3 && are_equal(buf, "inf"))
        return std::numeric_limits<Double>::infinity();

    constexpr std::size_t max_digits = 320;
    if (buf.size() > max_digits)
        return std::nullopt;

    char buffer[max_digits] = { '\0' };
    std::copy_n(buf.data(), buf.size(), std::begin(buffer));
    buffer[buf.size()] = '\0';

    auto posexp = buf.find("e-");
    if (posexp == std::string_view::npos) {
        posexp = buf.find("E-");
    }
    auto pos = buf.find('.');
    if (pos != std::string_view::npos) {
        int decimal = buf.substr(buf.find('.') + 1, (posexp == std::string_view::npos) ? posexp : (posexp - pos - 1)).size();
        if (posexp != std::string_view::npos) {
            decimal += atoi(buffer + posexp + 2);
        }
        //cout << "decimal:" << decimal << endl;
        if (baryonyx::precision < decimal) {
            baryonyx::precision = decimal;
        }
    } else if (posexp != std::string_view::npos) {
        int decimal = atoi(buffer + posexp + 2);
        //cout << "decimal:" << decimal << endl;
        if (baryonyx::precision < decimal) {
            baryonyx::precision = decimal;
        }
    }

    Double result = 0.0;
    if (auto read = std::sscanf(buffer, PrintFormatProb, &result); read > 0)
        return result;
    else
        return std::nullopt;
}

// constexpr
std::optional<real_token> read_real(const std::string_view buf_1,
    const std::string_view buf_2) noexcept
{
    std::optional<Double> value_opt;

    if (buf_1 == "-") {
        value_opt = read_real(buf_2);
        if (!value_opt)
            return real_token(-1.0, 1);
        else
            return real_token(-1.0 * *value_opt, 2);
    }

    if (buf_1 == "+") {
        value_opt = read_real(buf_2);
        if (!value_opt)
            return real_token(1.0, 1);
        else
            return real_token(*value_opt, 2);
    }

    value_opt = read_real(buf_1);
    if (!value_opt)
        return real_token(1.0, 0);

    return real_token(*value_opt, 1);
}

// constexpr
std::optional<std::string_view> read_name(const std::string_view buf) noexcept
{
    return std::find_if_not(std::begin(buf), std::end(buf), is_char_name) == std::end(buf)
        ? std::make_optional(buf)
        : std::nullopt;
}

// constexpr
std::optional<operator_token>
read_operator(const std::string_view buf_1,
    const std::string_view buf_2) noexcept
{
    if (buf_1 == "<")
        return operator_token{ baryonyx::operator_type::less, buf_2 == "=" ? 2 : 1 };

    if (buf_1 == ">")
        return operator_token{ baryonyx::operator_type::greater,
            buf_2 == "=" ? 2 : 1 };

    if (buf_1 == "=") {
        if (buf_2 == "<")
            return operator_token{ baryonyx::operator_type::less, 2 };

        if (buf_2 == "=")
            return operator_token{ baryonyx::operator_type::equal, 2 };

        if (buf_2 == ">")
            return operator_token{ baryonyx::operator_type::greater, 2 };

        return operator_token{ baryonyx::operator_type::equal, 1 };
    }

    return std::nullopt;
}

constexpr bool is_objective_function(const std::string_view buf) noexcept
{
    return !buf.empty() && !is_keyword(buf) && (buf == "+" || buf == "-" || is_number(buf[0]) || is_char_name(buf[0]));
}

std::optional<function_element_token>
read_quadratic_element(const std::string_view buf_1,
    const std::string_view buf_2,
    const std::string_view buf_3) noexcept
{
    auto value_opt = read_real(buf_1, buf_2);
    if (!value_opt)
        return std::nullopt;

    std::string_view to_read;
    switch (value_opt->read) {
    case 0:
        to_read = buf_1;
        break;
    case 1:
        to_read = buf_2;
        break;
    default:
        to_read = buf_3;
        break;
    }

    if (is_keyword(to_read) || !starts_with_name(to_read))
        return std::nullopt;

    auto name_opt = read_name(to_read);
    if (!name_opt.has_value())
        return std::nullopt;

    return function_element_token(value_opt->value, *name_opt,
        value_opt->read + 1);
}

raw_problem_status read_quadratic_element(stream_buffer& buf, problem_parser& p,
    Double sign_factor) noexcept
{
    if (buf.first() != "[")
        return baryonyx::file_format_error_tag::bad_objective_quadratic;

    buf.pop_front();

    while (!buf.first().empty() && buf.first() != "]") {
        function_element_token fct;

        if (auto fct_opt = read_quadratic_element(buf.first(), buf.second(), buf.third());
            !fct_opt)
            return baryonyx::file_format_error_tag::bad_objective_quadratic;
        else
            fct = *fct_opt;

        buf.pop_front(fct.read);

        if (buf.first() == "*") {
            auto name_opt = read_name(buf.second());
            if (!name_opt.has_value())
                return baryonyx::file_format_error_tag::bad_objective_quadratic;

            auto id1 = p.get_or_assign_variable(fct.name);
            auto id2 = p.get_or_assign_variable(*name_opt);

            if (id1 == -1 || id2 == -1)
                return baryonyx::file_format_error_tag::too_many_variables;

            auto it = std::find_if(p.m_problem.objective.qelements.begin(),
                p.m_problem.objective.qelements.end(),
                [id1, id2](const auto& elem) {
                    return (elem.variable_index_a == id1 && elem.variable_index_b == id2) || (elem.variable_index_a == id2 && elem.variable_index_b == id1);
                });
            if (it == p.m_problem.objective.qelements.end()) {
                auto& ret = p.m_problem.objective.qelements.emplace_back();
                ret.factor = fct.factor * sign_factor / 2.0;
                ret.variable_index_a = id1;
                ret.variable_index_b = id2;
            } else {
                it->factor += fct.factor * sign_factor / 2.0;
            }

            buf.pop_front(2);
        } else if (buf.first() == "^" || buf.first() == "^2") {
            if (buf.first() == "^" && buf.second() == "2")
                buf.pop_front(2);
            else
                buf.pop_front(1);

            auto id = p.get_or_assign_variable(fct.name);
            if (id == -1)
                return baryonyx::file_format_error_tag::too_many_variables;

            auto it = std::find_if(
                p.m_problem.objective.qelements.begin(),
                p.m_problem.objective.qelements.end(), [id](const auto& elem) {
                    return elem.variable_index_a == id && elem.variable_index_b == id;
                });

            if (it == p.m_problem.objective.qelements.end()) {
                auto& ret = p.m_problem.objective.qelements.emplace_back();
                ret.factor = fct.factor * sign_factor / 2.0;
                ret.variable_index_a = id;
                ret.variable_index_b = id;
            } else {
                it->factor += fct.factor * sign_factor / 2.0;
            }
        }
    }

    buf.pop_front();

    if (buf.first() == "/" && buf.second() == "2")
        buf.pop_front(2);
    else if (buf.first() == "/2")
        buf.pop_front(1);
    else
        return baryonyx::file_format_error_tag::bad_objective_quadratic;

    return baryonyx::file_format_error_tag::success;
}

// constexpr
std::optional<function_element_token>
read_function_element(const std::string_view buf_1,
    const std::string_view buf_2,
    const std::string_view buf_3) noexcept
{
    auto value_opt = read_real(buf_1, buf_2);
    if (!value_opt)
        return std::nullopt;

    std::string_view to_read;

    switch (value_opt->read) {
    case 0:
        to_read = buf_1;
        break;
    case 1:
        to_read = buf_2;
        break;
    default:
        to_read = buf_3;
        break;
    }

    if (!is_keyword(to_read) && starts_with_name(to_read)) {
        auto name_opt = read_name(to_read);
        if (!name_opt.has_value())
            return std::nullopt;

        return function_element_token(value_opt->value, *name_opt,
            value_opt->read + 1);
    } else {
        return function_element_token(value_opt->value, std::string_view(),
            value_opt->read);
    }
}

std::optional<sub_bound_token> read_left_bound_token(
    const std::string_view buf_1, const std::string_view buf_2,
    const std::string_view buf_3, const std::string_view buf_4) noexcept
{
    Double value = 0.0;
    Double negative = 1.0;

    if (buf_1 == "+" || buf_1 == "-") {
        if (buf_1 == "-")
            negative = -1.0;

        if (auto value_opt = read_real(buf_2); value_opt)
            value = *value_opt;
        else
            return std::nullopt;

        if (auto operator_token = read_operator(buf_3, buf_4); operator_token)
            return sub_bound_token(value * negative, operator_token->read + 2,
                operator_token->op);
        else
            return std::nullopt;
    } else {
        if (auto value_opt = read_real(buf_1); value_opt)
            value = *value_opt;
        else
            return std::nullopt;

        if (auto operator_token = read_operator(buf_2, buf_3); operator_token)
            return sub_bound_token(value * negative, operator_token->read + 1,
                operator_token->op);
        else
            return std::nullopt;
    }
}

std::optional<sub_bound_token> read_right_bound_token(
    const std::string_view buf_1, const std::string_view buf_2,
    const std::string_view buf_3, const std::string_view buf_4) noexcept
{
    Double negative = 1.0;
    baryonyx::operator_type op = baryonyx::operator_type::equal;

    if (auto operator_token = read_operator(buf_1, buf_2); operator_token) {
        op = operator_token->op;
        if (operator_token->read == 1) {
            if (buf_2 == "+" || buf_2 == "-") {
                if (buf_2 == "-")
                    negative = -1.0;

                if (auto value_opt = read_real(buf_3); value_opt)
                    return sub_bound_token(negative * *value_opt, 3, op);
                else
                    return std::nullopt;
            } else {
                if (auto value_opt = read_real(buf_2); value_opt)
                    return sub_bound_token(negative * *value_opt, 2, op);
                else
                    return std::nullopt;
            }
        } else {
            if (buf_3 == "+" || buf_3 == "-") {
                if (buf_3 == "-")
                    negative = -1.0;

                if (auto value_opt = read_real(buf_4); value_opt)
                    return sub_bound_token(negative * *value_opt, 4, op);
                else
                    return std::nullopt;
            } else {
                if (auto value_opt = read_real(buf_3); value_opt)
                    return sub_bound_token(negative * *value_opt, 3, op);
                else
                    return std::nullopt;
            }
        }
    } else
        return std::nullopt;
}

std::optional<bound_token>
read_bound(const stream_buffer::string_view_array& tokens) noexcept
{
    if (starts_with_number(tokens[0])) {
        auto left_opt = read_left_bound_token(tokens[0], tokens[1], tokens[2], tokens[3]);
        if (!left_opt)
            return std::nullopt;

        int read = left_opt->read;

        auto name_opt = read_name(tokens[read]);
        if (!name_opt)
            return std::nullopt;

        ++read;

        if (!starts_with_operator(tokens[read]))
            return bound_token(left_opt->value,
                std::numeric_limits<Double>::infinity(), *name_opt,
                read + 1);

        auto right_opt = read_right_bound_token(tokens[read], tokens[read + 1],
            tokens[read + 2], tokens[read + 3]);
        if (!right_opt)
            return std::nullopt;

        if (left_opt->value > right_opt->value)
            return std::nullopt;

        return bound_token(left_opt->value, right_opt->value, *name_opt,
            read + right_opt->read);
    }

    if (starts_with_name(tokens[0])) {
        auto name_opt = read_name(tokens[0]);
        if (!name_opt)
            return std::nullopt;

        if (starts_with_operator(tokens[1])) {
            auto right_opt = read_right_bound_token(tokens[1], tokens[2], tokens[3], tokens[4]);
            if (!right_opt)
                return std::nullopt;

            if (tokens[1] == "=") {
                return bound_token(right_opt->value, right_opt->value, *name_opt, 1 + right_opt->read);
            } else if (tokens[1] == ">") {
                return bound_token(right_opt->value, std::numeric_limits<Double>::infinity(), *name_opt, 1 + right_opt->read);
            } else if (tokens[1] == "<") {
                return bound_token(0.0, right_opt->value, *name_opt, 1 + right_opt->read);
            } else {
                return std::nullopt;
            }
        }

        return bound_token(-std::numeric_limits<Double>::infinity(),
            +std::numeric_limits<Double>::infinity(), *name_opt, 1);
    }

    return std::nullopt;
}

std::optional<baryonyx::objective_function_type>
read_objective_type(const std::string_view token) noexcept
{
    if (are_equal(token, "maximize") || are_equal(token, "maximum") || are_equal(token, "max"))
        return baryonyx::objective_function_type::maximize;

    if (are_equal(token, "minimize") || are_equal(token, "minimum") || are_equal(token, "min"))
        return baryonyx::objective_function_type::minimize;

    return std::nullopt;
}

constexpr int try_read_objective_label(const std::string_view first,
    const std::string_view second) noexcept
{
    if (is_keyword(first))
        return 0;

    if (are_equal(second, ":"))
        return 2;

    return 0;
}

constexpr int read_subject_to(const std::string_view buf_1,
    const std::string_view buf_2,
    const std::string_view buf_3) noexcept
{
    if (are_equal(buf_1, "st") || are_equal(buf_1, "st.") || are_equal(buf_1, "s.t") || are_equal(buf_1, "s.t."))
        return are_equal(buf_2, ":") ? 2 : 1;

    if (are_equal(buf_1, "subject"))
        if (are_equal(buf_2, "to"))
            return are_equal(buf_3, ":") ? 3 : 2;

    if (are_equal(buf_1, "sush"))
        if (are_equal(buf_2, "that"))
            return are_equal(buf_3, ":") ? 3 : 2;

    return 0;
}

constexpr int read_bounds(const std::string_view buf_1,
    const std::string_view buf_2) noexcept
{
    if (are_equal(buf_1, "bounds") || are_equal(buf_1, "bound"))
        return are_equal(buf_2, ":") ? 2 : 1;

    return 0;
}

constexpr int read_binary(const std::string_view buf_1,
    const std::string_view buf_2) noexcept
{
    if (are_equal(buf_1, "binary") || are_equal(buf_1, "binaries") || are_equal(buf_1, "bin"))
        return are_equal(buf_2, ":") ? 2 : 1;

    return 0;
}

constexpr int read_general(const std::string_view buf_1,
    const std::string_view buf_2) noexcept
{
    if (are_equal(buf_1, "general") || are_equal(buf_1, "generals") || are_equal(buf_1, "gen"))
        return are_equal(buf_2, ":") ? 2 : 1;

    return 0;
}

constexpr int read_end(const std::string_view buf_1,
    const std::string_view buf_2) noexcept
{
    if (are_equal(buf_1, "end"))
        return are_equal(buf_2, ":") ? 2 : 1;

    return 0;
}

bool starts_with_quadratic(const std::string_view buf_1,
    const std::string_view buf_2) noexcept
{
    return buf_1 == "[" || (buf_1 == "+" && buf_2 == "[") || (buf_1 == "-" && buf_2 == "[");
}

raw_problem_status parse(stream_buffer& buf, problem_parser& p) noexcept
{
    baryonyx::precision = 0;
    int obj_precision = 0;

    if (auto obj_type = read_objective_type(buf.first()); !obj_type)
        return baryonyx::file_format_error_tag::bad_objective_function_type;
    else
        p.m_problem.type = *obj_type;

    buf.pop_front();

    if (auto read = try_read_objective_label(buf.first(), buf.second()); read)
        buf.pop_front(read);

    while (!is_keyword(buf.first())) {
        if (starts_with_quadratic(buf.first(), buf.second())) {
            Double factor = 1.0;

            if (buf.first() == "-") {
                factor = -1.0;
                buf.pop_front();
            } else if (buf.first() == "+") {
                buf.pop_front();
            }

            if (auto ret = read_quadratic_element(buf, p, factor); !ret)
                return ret;
            continue;
        }

        if (auto elem = read_function_element(buf.first(), buf.second(), buf.third());
            elem) {
            if (!p.append_to_objective(*elem))
                return baryonyx::file_format_error_tag::too_many_variables;

            buf.pop_front(elem->read);
            continue;
        }

        return baryonyx::file_format_error_tag::bad_objective;
    }

    obj_precision = baryonyx::precision;

    if (auto read = read_subject_to(buf.first(), buf.second(), buf.third());
        read) {
        buf.pop_front(read);
        int constraint_next_id = 0;

        while (!buf.first().empty() && !is_keyword(buf.first())) {
            std::string_view label;
            Double value;
            std::vector<baryonyx::function_element> elements;
            baryonyx::precision = 0;

            if (starts_with_name(buf.first()) && buf.second() == ":") {
                label = p.m_problem.strings->append(buf.first());
                buf.pop_front(2);
            }

            auto elem = read_function_element(buf.first(), buf.second(), buf.third());
            if (!elem)
                return baryonyx::file_format_error_tag::bad_constraint;

            {
                auto id = p.get_or_assign_variable(elem->name);
                if (id == -1)
                    return baryonyx::file_format_error_tag::too_many_variables;

                auto it = std::find_if(elements.begin(), elements.end(),
                    [id](const auto& func_elem) {
                        return func_elem.variable_index == id;
                    });

                if (it == elements.end())
                    elements.emplace_back(elem->factor, id);
                else
                    it->factor += elem->factor;
            }

            buf.pop_front(elem->read);

            while (!buf.first().empty() && !starts_with_operator(buf.first())) {
                elem = read_function_element(buf.first(), buf.second(), buf.third());
                if (!elem)
                    return baryonyx::file_format_error_tag::bad_constraint;

                {
                    auto id = p.get_or_assign_variable(elem->name);
                    if (id == -1)
                        return baryonyx::file_format_error_tag::too_many_variables;

                    auto it = std::find_if(elements.begin(), elements.end(),
                        [id](const auto& func_elem) {
                            return func_elem.variable_index == id;
                        });

                    if (it == elements.end())
                        elements.emplace_back(elem->factor, id);
                    else
                        it->factor += elem->factor;
                }

                if (elements.back().variable_index == -1)
                    return baryonyx::file_format_error_tag::too_many_variables;

                buf.pop_front(elem->read);
            }

            auto operator_opt = read_operator(buf.first(), buf.second());
            if (!operator_opt)
                return baryonyx::file_format_error_tag::bad_constraint;
            buf.pop_front(operator_opt->read);

            auto value_opt = read_real(buf.first(), buf.second());
            if (!value_opt)
                return baryonyx::file_format_error_tag::bad_constraint;
            value = value_opt->value;
            buf.pop_front(value_opt->read);

            switch (operator_opt->op) {
            case baryonyx::operator_type::equal:
                if (elements.size() == 1) {
                    int val = static_cast<int>(baryonyx::Floor(value / elements[0].factor));
                    assert(val >= p.m_problem.vars.values[elements[0].variable_index].min);
                    assert(val <= p.m_problem.vars.values[elements[0].variable_index].max);
                    p.m_problem.vars.values[elements[0].variable_index].min = val;
                    p.m_problem.vars.values[elements[0].variable_index].max = val;
                } else {
                    p.m_problem.equal_constraints.emplace_back(label, std::move(elements), value, constraint_next_id++, baryonyx::precision);
                }
                break;
            case baryonyx::operator_type::greater:
                if (elements.size() == 1) {
                    int val = static_cast<int>(baryonyx::Ceil(value / elements[0].factor));
                    assert(val <= p.m_problem.vars.values[elements[0].variable_index].max);
                    p.m_problem.vars.values[elements[0].variable_index].min = max(p.m_problem.vars.values[elements[0].variable_index].min, val);
                } else {
                    p.m_problem.greater_constraints.emplace_back(label, std::move(elements), value, constraint_next_id++, baryonyx::precision);
                }
                break;
            case baryonyx::operator_type::less:
                if (elements.size() == 1) {
                    int val = static_cast<int>(baryonyx::Floor(value / elements[0].factor));
                    assert(val >= p.m_problem.vars.values[elements[0].variable_index].min);
                    p.m_problem.vars.values[elements[0].variable_index].max = min(p.m_problem.vars.values[elements[0].variable_index].max, val);
                } else {
                    p.m_problem.less_constraints.emplace_back(label, std::move(elements), value, constraint_next_id++, baryonyx::precision);
                }
                break;
            }
        }

        if (std::cmp_greater(p.m_problem.equal_constraints.size() + p.m_problem.greater_constraints.size() + p.m_problem.less_constraints.size(),
                INT_MAX))
            return baryonyx::file_format_error_tag::too_many_constraints;
    }

    if (auto read = read_bounds(buf.first(), buf.second()); read) {
        buf.pop_front(read);

        while (!is_keyword(buf.first())) {
            auto bound_opt = read_bound(buf.array());
            if (!bound_opt)
                return baryonyx::file_format_error_tag::bad_bound;

            if (!p.set_bound_variable(*bound_opt))
                return baryonyx::file_format_error_tag::bad_bound;

            buf.pop_front(bound_opt->read);
        }
    }

    if (auto read = read_binary(buf.first(), buf.second()); read) {
        buf.pop_front(read);

        while (!is_keyword(buf.first()))
            if (!p.set_boolean_variable(buf.first()))
                return baryonyx::file_format_error_tag::bad_binary;
            else
                buf.pop_front();
    }

    if (auto read = read_general(buf.first(), buf.second()); read) {
        buf.pop_front(read);

        while (!is_keyword(buf.first()))
            if (!p.set_integer_variable(buf.first()))
                return baryonyx::file_format_error_tag::bad_general;
            else
                buf.pop_front();
    }

    baryonyx::precision = obj_precision;

    if (auto read = read_end(buf.first(), buf.second()); read) {
        buf.pop_front(read);

        return buf.first().empty()
            ? baryonyx::file_format_error_tag::success
            : baryonyx::file_format_error_tag::bad_end_of_file;
    }

    return baryonyx::file_format_error_tag::bad_end;
}

namespace baryonyx {

raw_problem make_problem(std::istream& is) noexcept
{
    stream_buffer sb(is);
    problem_parser p;

    if (auto ret = parse(sb, p); ret)
        return std::move(p.m_problem);
    else
        return raw_problem(ret.tag);
}

raw_problem make_problem(std::string_view file_path) noexcept
{
    const std::string copy{ file_path };
    std::ifstream ifs{ copy };

    return make_problem(ifs);
}

} // namespace baryonyx
