#!/usr/bin/env python3
"""
Entry point script for decoy generation package.

This script can be run from anywhere and properly handles the package imports.
"""

import sys
import os
from pathlib import Path

# Add the package directory to Python path
package_dir = Path(__file__).parent
sys.path.insert(0, str(package_dir))

# Now we can import and run the main function
if __name__ == "__main__":
    from main import main
    main()