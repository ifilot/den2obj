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

#ifndef _CLI_FORMAT
#define _CLI_FORMAT

#include <cstddef>
#include <iostream>
#include <string>

namespace CLIFormat {

    enum class Color {
        DEFAULT,
        BOLD,
        DIM,
        CYAN,
        GREEN,
        YELLOW
    };

    /**
     * @brief      Format a byte count as kibibytes.
     *
     * @param[in]  bytes  Number of bytes
     *
     * @return     Formatted size string
     */
    std::string format_kb(size_t bytes);

    /**
     * @brief      Format an elapsed time in seconds.
     *
     * @param[in]  seconds  Elapsed time in seconds
     *
     * @return     Formatted elapsed time string
     */
    std::string format_seconds(double seconds);

    /**
     * @brief      Check whether standard output is attached to a terminal.
     *
     * @return     Whether standard output is a terminal
     */
    bool output_is_terminal();

    /**
     * @brief      Check whether colored output should be emitted.
     *
     * @return     Whether ANSI color output is enabled
     */
    bool color_enabled();

    /**
     * @brief      Wrap a string in ANSI color/style codes when enabled.
     *
     * @param[in]  value  Text to decorate
     * @param[in]  color  Color or style to apply
     *
     * @return     Decorated text, or the original text when color is disabled
     */
    std::string colorize(const std::string& value, Color color);

    /**
     * @brief      Print the program banner.
     *
     * @param[in]  version          Program version
     * @param[in]  git_hash         Git commit hash
     * @param[in]  build_date       Compilation date
     * @param[in]  build_time       Compilation time
     * @param[in]  openvdb_enabled  Whether OpenVDB support is enabled
     */
    void print_banner(const std::string& version,
                      const std::string& git_hash,
                      const std::string& build_date,
                      const std::string& build_time,
                      bool openvdb_enabled);

    /**
     * @brief      Print a section heading.
     *
     * @param[in]  title  Section title
     */
    void print_section(const std::string& title);

    /**
     * @brief      Print a left-aligned key-value line.
     *
     * @param[in]  key    Field name
     * @param[in]  value  Field value
     * @param      out    Output stream
     */
    void print_kv(const std::string& key,
                  const std::string& value,
                  std::ostream& out = std::cout);

    /**
     * @brief      Print the final completion summary.
     *
     * @param[in]  elapsed  Formatted elapsed time
     */
    void print_done_summary(const std::string& elapsed);

    /**
     * @brief      Lightweight terminal-aware progress bar.
     */
    class ProgressBar {
    public:
        /**
         * @brief      Construct a progress bar.
         *
         * @param[in]  label  Progress label
         * @param[in]  total  Total number of ticks
         */
        ProgressBar(const std::string& label, size_t total);

        /**
         * @brief      Destroy the progress bar and finish its line.
         */
        ~ProgressBar();

        /**
         * @brief      Advance the progress bar by one tick.
         */
        void tick();

    private:
        /**
         * @brief      Render the progress bar.
         *
         * @param[in]  force  Whether to render even if percentage did not change
         */
        void render(bool force);

        std::string label;
        size_t total;
        size_t current;
        size_t last_percent;
        bool terminal;
    };

}

#endif
