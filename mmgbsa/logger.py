#!/usr/bin/env python3
"""
Custom Logger for OpenGBSA
Provides professional, rich-formatted output for CLI tools
"""
import sys
import logging
from datetime import datetime

class ToolLogger:
    """
    Color-coded logger for CLI output
    """
    
    # ANSI Color Codes
    BLUE = '\033[94m'
    CYAN = '\033[96m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    ENDC = '\033[0m'
    
    # Icons
    ICON_INFO = "ℹ️ "
    ICON_SUCCESS = "✅"
    ICON_WARNING = "⚠️ "
    ICON_ERROR = "❌"
    ICON_PROCESS = "⚙️ "
    ICON_TIME = "⏱️ "
    ICON_FILE = "📁"
    ICON_ROCKET = "🚀"
    
    def __init__(self, name="OpenGBSA", verbose=1):
        self.name = name
        self.verbose = verbose
        
    def header(self, message):
        """Print a major section header"""
        print(f"\n{self.BOLD}{self.BLUE}" + "="*80 + f"{self.ENDC}")
        print(f"{self.BOLD}{self.BLUE}{self.ICON_ROCKET} {message.center(76)} {self.ENDC}")
        print(f"{self.BOLD}{self.BLUE}" + "="*80 + f"{self.ENDC}\n")
        
    def section(self, message):
        """Print a subsection header"""
        print(f"\n{self.BOLD}{self.CYAN}--- {message} ---{self.ENDC}")

    def info(self, message):
        """Print informational message"""
        if self.verbose >= 1:
            print(f"{self.BLUE}{self.ICON_INFO} {message}{self.ENDC}")
            
    def success(self, message):
        """Print success message"""
        print(f"{self.GREEN}{self.ICON_SUCCESS} {message}{self.ENDC}")
        
    def warning(self, message):
        """Print warning message"""
        print(f"{self.YELLOW}{self.ICON_WARNING} {message}{self.ENDC}")
        
    def error(self, message):
        """Print error message"""
        print(f"{self.RED}{self.BOLD}{self.ICON_ERROR} {message}{self.ENDC}", file=sys.stderr)
        
    def process(self, message):
        """Print processing status"""
        if self.verbose >= 1:
            print(f"{self.CYAN}{self.ICON_PROCESS} {message}{self.ENDC}")
            
    def result(self, label, value, unit=""):
        """Print a result value"""
        print(f"  • {label}: {self.BOLD}{value}{self.ENDC} {unit}")
        
    def debug(self, message):
        """Print debug message"""
        if self.verbose >= 2:
            print(f"{self.YELLOW}[DEBUG] {message}{self.ENDC}")

logger = ToolLogger()
