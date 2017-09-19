"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Unicycler

This module contains a Unicycler class for writing output to both stdout and a log file.

This file is part of Unicycler. Unicycler is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Unicycler is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Unicycler. If
not, see <http://www.gnu.org/licenses/>.
"""

import sys
import datetime
import re
import shutil
import textwrap
import subprocess


class Log(object):

    def __init__(self):
        try:
            self.colours = int(subprocess.check_output(['tput', 'colors']).decode().strip())
        except (ValueError, subprocess.CalledProcessError, FileNotFoundError, AttributeError):
            self.colours = 1


# This is the one and only instance of the Log class.
logger = Log()


def log(text, end='\n'):
    text_no_formatting = remove_formatting(text)

    # The text is printed to the screen with ANSI formatting, if supported. If there are only 8
    # colours available, then remove the 'dim' format which doesn't work.
    if logger.colours <= 1:
        text = text_no_formatting
    elif logger.colours <= 8:
        text = remove_dim_formatting(text)
    print(text, file=sys.stderr, end=end, flush=True)


def log_section_header(message, single_newline=False):
    if single_newline:
        log('')
    else:
        log('\n')

    time = get_timestamp()
    time_str = '(' + time + ')'
    if logger.colours > 8:
        time_str = dim(time_str)
    log(bold_yellow_underline(message) + ' ' + time_str)


def log_explanation(text, extra_empty_lines_after=1, indent_size=4):
    text = ' ' * indent_size + text
    terminal_width = shutil.get_terminal_size().columns
    for line in textwrap.wrap(text, width=terminal_width - 1):
        if logger.colours > 8:
            formatted_text = dim(line)
        else:
            formatted_text = line
        log(formatted_text)

    for _ in range(extra_empty_lines_after):
        log('')


def get_timestamp():
    return '{:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now())


END_FORMATTING = '\033[0m'
BOLD = '\033[1m'
UNDERLINE = '\033[4m'
YELLOW = '\033[93m'
DIM = '\033[2m'


def bold_yellow_underline(text):
    return YELLOW + BOLD + UNDERLINE + text + END_FORMATTING


def dim(text):
    return DIM + text + END_FORMATTING


def remove_formatting(text):
    return re.sub('\033.*?m', '', text)


def remove_dim_formatting(text):
    return re.sub('\033\[2m', '', text)
