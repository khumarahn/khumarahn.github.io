#pragma once

#include <stdexcept>
#include <string_view>
#include <source_location>
#include <format>

inline void verify(bool condition, 
        std::string_view message = "error here", 
        const std::source_location loc = std::source_location::current()) 
{
    if (!condition) {
        // std::format constructs a clean std::string, which perfectly 
        // satisfies std::runtime_error's constructor.
        throw std::runtime_error(
                std::format("{}:{}: {}\n", loc.file_name(), loc.line(), message)
                );
    }
}
