/**************************************************************************
 *                                                                        *
 *   Author: Ivo Filot <i.a.w.filot@tue.nl>                               *
 *                                                                        *
 *   DEN2OBJ is free software:                                            *
 *   you can redistribute it and/or modify it under the terms of the      *
 *   GNU General Public License as published by the Free Software         *
 *   Foundation, either version 3 of the License, or (at your option)     *
 *   any later version.                                                   *
 *                                                                        *
 *   DEN2OBJ is distributed in the hope that it will be useful,           *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty          *
 *   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.              *
 *   See the GNU General Public License for more details.                 *
 *                                                                        *
 *   You should have received a copy of the GNU General Public License    *
 *   along with this program.  If not, see http://www.gnu.org/licenses/.  *
 *                                                                        *
 **************************************************************************/

#include "cli_format.h"

#include <algorithm>
#include <cstdlib>
#include <iomanip>
#include <sstream>
#include <unistd.h>

namespace {
    /**
     * @brief      Pad a string on the right with spaces.
     *
     * @param[in]  value  Value to pad
     * @param[in]  width  Minimum output width
     *
     * @return     Padded string
     */
    std::string pad_right(const std::string& value, size_t width) {
        if(value.size() >= width) {
            return value;
        }

        return value + std::string(width - value.size(), ' ');
    }

    /**
     * @brief      Get the ANSI escape sequence for a color or style.
     *
     * @param[in]  color  Color or style
     *
     * @return     ANSI escape sequence
     */
    std::string ansi_code(CLIFormat::Color color) {
        switch(color) {
            case CLIFormat::Color::BOLD:
                return "\033[1m";
            case CLIFormat::Color::DIM:
                return "\033[2m";
            case CLIFormat::Color::CYAN:
                return "\033[36m";
            case CLIFormat::Color::GREEN:
                return "\033[32m";
            case CLIFormat::Color::YELLOW:
                return "\033[33m";
            case CLIFormat::Color::DEFAULT:
            default:
                return "\033[0m";
        }
    }
}

/**
 * @brief      Format a byte count as kibibytes.
 *
 * @param[in]  bytes  Number of bytes
 *
 * @return     Formatted size string
 */
std::string CLIFormat::format_kb(size_t bytes) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(1) << (static_cast<double>(bytes) / 1024.0) << " kb";
    return oss.str();
}

/**
 * @brief      Format an elapsed time in seconds.
 *
 * @param[in]  seconds  Elapsed time in seconds
 *
 * @return     Formatted elapsed time string
 */
std::string CLIFormat::format_seconds(double seconds) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(2) << seconds << " s";
    return oss.str();
}

/**
 * @brief      Check whether standard output is attached to a terminal.
 *
 * @return     Whether standard output is a terminal
 */
bool CLIFormat::output_is_terminal() {
    return isatty(STDOUT_FILENO) != 0;
}

/**
 * @brief      Check whether colored output should be emitted.
 *
 * @return     Whether ANSI color output is enabled
 */
bool CLIFormat::color_enabled() {
    return CLIFormat::output_is_terminal() && std::getenv("NO_COLOR") == nullptr;
}

/**
 * @brief      Wrap a string in ANSI color/style codes when enabled.
 *
 * @param[in]  value  Text to decorate
 * @param[in]  color  Color or style to apply
 *
 * @return     Decorated text, or the original text when color is disabled
 */
std::string CLIFormat::colorize(const std::string& value, Color color) {
    if(!CLIFormat::color_enabled()) {
        return value;
    }

    return ansi_code(color) + value + ansi_code(Color::DEFAULT);
}

/**
 * @brief      Print the program banner.
 *
 * @param[in]  version          Program version
 * @param[in]  git_hash         Git commit hash
 * @param[in]  build_date       Compilation date
 * @param[in]  build_time       Compilation time
 * @param[in]  openvdb_enabled  Whether OpenVDB support is enabled
 */
void CLIFormat::print_banner(const std::string& version,
                             const std::string& git_hash,
                             const std::string& build_date,
                             const std::string& build_time,
                             bool openvdb_enabled) {
    static const size_t width = 48;
    const std::string rule = "+" + std::string(width, '-') + "+";

    std::cout << colorize(rule, Color::CYAN) << std::endl;
    std::cout << colorize("| " + pad_right("DEN2OBJ " + version, width - 2) + " |", Color::BOLD) << std::endl;
    std::cout << "| " << pad_right("Density fields -> meshes and volumes", width - 2) << " |" << std::endl;
    std::cout << colorize(rule, Color::CYAN) << std::endl;
    print_kv("Author", "Ivo Filot <i.a.w.filot@tue.nl>");
    print_kv("Website", "https://den2obj.imc-tue.nl");
    print_kv("GitHub", "https://github.com/ifilot/den2obj");
    print_kv("Build", build_date + " " + build_time);
    print_kv("Git", git_hash);
    if(openvdb_enabled) {
        print_kv("Modules", "OpenVDB");
    }
    std::cout << std::endl;
}

/**
 * @brief      Print a section heading.
 *
 * @param[in]  title  Section title
 */
void CLIFormat::print_section(const std::string& title) {
    std::cout << std::endl << colorize("== " + title + " ==", Color::CYAN) << std::endl;
}

/**
 * @brief      Print a left-aligned key-value line.
 *
 * @param[in]  key    Field name
 * @param[in]  value  Field value
 * @param      out    Output stream
 */
void CLIFormat::print_kv(const std::string& key,
                         const std::string& value,
                         std::ostream& out) {
    out << colorize(pad_right(key + ":", 14), Color::DIM) << value << std::endl;
}

/**
 * @brief      Print the final completion summary.
 *
 * @param[in]  elapsed  Formatted elapsed time
 */
void CLIFormat::print_done_summary(const std::string& elapsed) {
    print_section("Done");
    print_kv("Elapsed", elapsed);
    std::cout << std::endl;
}

/**
 * @brief      Construct a progress bar.
 *
 * @param[in]  label  Progress label
 * @param[in]  total  Total number of ticks
 */
CLIFormat::ProgressBar::ProgressBar(const std::string& label, size_t total) :
    label(label),
    total(total),
    current(0),
    last_percent(static_cast<size_t>(-1)),
    terminal(CLIFormat::output_is_terminal()) {
    if(this->terminal) {
        this->render(true);
    }
}

/**
 * @brief      Destroy the progress bar and finish its line.
 */
CLIFormat::ProgressBar::~ProgressBar() {
    if(this->current < this->total) {
        this->current = this->total;
    }
    if(!this->terminal || this->last_percent != 100) {
        this->render(true);
    }
    std::cout << std::endl;
}

/**
 * @brief      Advance the progress bar by one tick.
 */
void CLIFormat::ProgressBar::tick() {
    if(this->current < this->total) {
        this->current++;
    }
    this->render(false);
}

/**
 * @brief      Render the progress bar.
 *
 * @param[in]  force  Whether to render even if percentage did not change
 */
void CLIFormat::ProgressBar::render(bool force) {
    static const size_t width = 40;
    const size_t percent = this->total == 0 ? 100 : (this->current * 100) / this->total;

    if(!force && percent == this->last_percent) {
        return;
    }

    if(!this->terminal && (percent != 100 || !force)) {
        return;
    }

    const size_t filled = std::min(width, (percent * width) / 100);
    const std::string bar = CLIFormat::colorize(std::string(filled, '#'), Color::GREEN) +
                            std::string(width - filled, '-');

    if(this->terminal) {
        std::cout << "\r";
    }

    std::cout << this->label << " [" << bar << "] "
              << std::right << std::setw(3) << percent << "%" << std::flush;

    this->last_percent = percent;
}
