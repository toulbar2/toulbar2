/* Copyright (C) 2016-2023 INRAE @ Gauthier Quesnel
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

#ifndef ORG_VLEPROJECT_BARYONYX_SOLVER_CORE
#define ORG_VLEPROJECT_BARYONYX_SOLVER_CORE

#include <array>
#include <forward_list>
#include <functional>
#include <limits>
#include <memory>
#include <numeric>
#include <string>
#include <string_view>
#include <vector>
#include <algorithm>

#include <cstdint>
#include <cmath>

#include <type_traits>

typedef long double Double;

namespace baryonyx {

/** @c index is used as accessors for all array. */
using index = std::int32_t;

/** @c value is used as value for variable. */
using var_value = std::int_least8_t;

Double Floor(Double v);
Double Ceil(Double v);

struct string_buffer {
    constexpr static std::size_t string_buffer_node_length = 1024 * 1024;

    using value_type = std::array<char, string_buffer_node_length>;
    using container_type = std::forward_list<value_type>;

    string_buffer() noexcept = default;
    ~string_buffer() noexcept = default;

    string_buffer(const string_buffer&) = delete;
    string_buffer(string_buffer&&) = delete;
    string_buffer& operator=(const string_buffer&) = delete;
    string_buffer& operator=(string_buffer&&) = delete;

    std::string_view append(std::string_view str)
    {
        if (m_container.empty() || str.size() + m_position > string_buffer_node_length)
            do_alloc();

        std::size_t position = m_position;
        m_position += str.size();

        char* buffer = m_container.front().data() + position;

        std::copy_n(str.data(), str.size(), buffer);

        return std::string_view(buffer, str.size());
    }

private:
    void do_alloc()
    {
        m_container.emplace_front();
        m_position = 0;
    }

    container_type m_container;
    std::size_t m_position = { 0 };
};

using string_buffer_ptr = std::shared_ptr<string_buffer>;

enum class file_format_error_tag {
    success,
    file_not_found,
    bad_end_of_file,
    bad_general,
    bad_binary,
    bad_objective_function_type,
    bad_objective,
    bad_objective_quadratic,
    bad_bound,
    bad_end,
    bad_constraint,
    too_many_variables,
    too_many_constraints,
    bad_name, // Mainly use when read result file and report unknown variable.
    empty_context
};

enum class variable_type { real,
    binary,
    general };

enum class objective_function_type { maximize,
    minimize };

enum class operator_type {
    equal,
    greater,
    less,
};

struct variable_value {
    variable_value() = default;

    bool touched{ false }; // true if min value is explicitly specified in the problem file (either a unary constraint or a left bound or a binary variable)
    int min{ std::numeric_limits<int>::min() };
    int max{ std::numeric_limits<int>::max() };
    variable_type type{ variable_type::real };
};

struct variables {
    std::vector<std::string_view> names;
    std::vector<variable_value> values;
};

struct function_element {
    function_element() = default;

    function_element(Double factor_, index variable_index_) noexcept
        : factor(factor_)
        , variable_index(variable_index_)
    {
    }

    Double factor = { 0.0 };
    index variable_index{ -1 };
};

struct objective_function_element {
    objective_function_element(Double factor_, index variable_index_) noexcept
        : factor(factor_)
        , variable_index(variable_index_)
    {
    }

    Double factor = { 0.0 };
    index variable_index{ -1 };
};

struct objective_quadratic_element {
    objective_quadratic_element() noexcept = default;

    constexpr objective_quadratic_element(Double factor_, index variable_index_a_,
        index variable_index_b_) noexcept
        : factor(factor_)
        , variable_index_a(variable_index_a_)
        , variable_index_b(variable_index_b_)
    {
    }

    Double factor = { 0.0 };
    index variable_index_a{ -1 };
    index variable_index_b{ -1 };
};

struct constraint {
    constraint() = default;

    constraint(std::string_view label_, std::vector<function_element>&& elements_,
        Double value_, int id_, int precision_)
        : label(label_)
        , elements(elements_)
        , value(value_)
        , id(id_)
        , precision(precision_)
    {
    }

    std::string_view label;
    std::vector<function_element> elements;
    Double value = { 0.0 };
    int id;
    int precision = { 0 };
};

struct objective_function {
    std::vector<objective_function_element> elements;
    std::vector<objective_quadratic_element> qelements;
    Double value = { 0.0 };
};

struct affected_variables {
    void push_back(std::string_view name, bool value)
    {
        names.emplace_back(name);
        values.emplace_back(value);
    }

    std::vector<std::string_view> names;
    std::vector<var_value> values;
};

struct raw_problem {
    raw_problem() = default;
    raw_problem(file_format_error_tag status_)
        : status(status_)
    {
    }

    string_buffer_ptr strings;

    objective_function objective;

    std::vector<constraint> equal_constraints;
    std::vector<constraint> greater_constraints;
    std::vector<constraint> less_constraints;

    constexpr operator bool() const noexcept
    {
        return status == file_format_error_tag::success;
    }

    variables vars;

    objective_function_type type = { objective_function_type::maximize };
    file_format_error_tag status = { file_format_error_tag::success };
};

raw_problem make_problem(std::istream& is) noexcept;

raw_problem make_problem(std::string_view file_path) noexcept;

extern int precision;

} // namespace baryonyx

std::ostream& operator<<(std::ostream& os, const baryonyx::file_format_error_tag& obj);

#endif
