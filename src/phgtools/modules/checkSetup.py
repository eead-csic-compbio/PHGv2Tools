#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
checkSetup.py
Check if all dependencies for PHGv2Tools are properly installed.
"""

import subprocess
import sys
import argparse


def check_phgtools():
    """Check if phgtools CLI is callable."""
    try:
        result = subprocess.run(['phgtools', '--version'], capture_output=True, text=True, check=True)
        print("✓ phgtools:", result.stdout.strip())
        return True
    except subprocess.CalledProcessError:
        print("✗ phgtools: not callable. Please check your PATH.")
        return False
    except FileNotFoundError:
        print("✗ phgtools: not installed or not found in PATH.")
        return False


def check_samtools():
    """Check if samtools is callable (required by fastaFromKey module)."""
    try:
        result = subprocess.run(['samtools', '--version'], capture_output=True, text=True, check=True)
        version_line = result.stdout.split('\n')[0]
        print(f"✓ samtools: {version_line}")
        return True
    except subprocess.CalledProcessError:
        print("✗ samtools: not callable. Please check your PATH.")
        return False
    except FileNotFoundError:
        print("✗ samtools: not installed or not found in PATH.")
        return False


def check_python():
    """Check Python version."""
    version = f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
    print(f"✓ python: {version}")
    return True


def check_python_packages():
    """Check required Python packages."""
    packages = {
        'matplotlib': 'matplotlib',
        'pandas': 'pandas',
        'scipy': 'scipy',
        'numpy': 'numpy',
        'tqdm': 'tqdm',
        'openpyxl': 'openpyxl',
        'rich': 'rich',
        'seaborn': 'seaborn',
    }
    
    all_ok = True
    print("\nPython packages:")
    for display_name, import_name in packages.items():
        try:
            module = __import__(import_name)
            version = getattr(module, '__version__', 'unknown')
            print(f"  ✓ {display_name}: {version}")
        except ImportError:
            print(f"  ✗ {display_name}: not installed")
            all_ok = False
    
    return all_ok


def main(args=None):
    """Main function to check all dependencies."""
    if args is None:
        args = sys.argv[1:]

    parser = argparse.ArgumentParser(
        prog="phgtools check-setup",
        description="Check if all dependencies for PHGv2Tools are properly installed.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="\nThis command checks:\n"
               "  - Command-line tools: phgtools, samtools\n"
               "  - Python version\n"
               "  - Python packages: matplotlib, pandas, scipy, numpy, tqdm, openpyxl, rich, seaborn\n"
               "\nExample usage:\n"
               "  phgtools check-setup"
    )
    
    parser.add_argument("-v", "--verbose", action="store_true", 
                        help="Print additional information")

    parsed_args = parser.parse_args(args)

    print("=" * 50)
    print("PHGv2Tools - Dependency Check")
    print("=" * 50)
    
    print("\nCommand-line tools:")
    check_phgtools()
    check_samtools()
    check_python()
    check_python_packages()
    
    print("\n" + "=" * 50)
    print("Check complete!")
    print("=" * 50)


if __name__ == "__main__":
    main()
