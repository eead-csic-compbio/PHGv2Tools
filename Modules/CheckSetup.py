import subprocess

def check_phgtools():
    try:
        # Check if phgtools is callable
        print("checking phgtools")
        result = subprocess.run(['phgtools', '--version'], capture_output=True, text=True, check=True)
        print("phgtools version:")
        print(result.stdout)
    except subprocess.CalledProcessError:
        print("phgtools is not callable from everywhere. Please check your PATH.")
    except FileNotFoundError:
        print("phgtools is not installed or not found in the PATH.")

def check_phgv2():
    try:
        # Check if phg is callable
        result = subprocess.run(['phg', '--version'], capture_output=True, text=True, check=True)
        print("checking phgv2")
        print(result.stdout)
    except subprocess.CalledProcessError:
        print("phg is not callable from everywhere. Please check your PATH.")
    except FileNotFoundError:
        print("phg is not installed or not found in the PATH.")

def check_samtools():
    try:
        # Check if samtools is callable
        result = subprocess.run(['samtools', '--version'], capture_output=True, text=True, check=True)
        print("checking samtools")
        print(result.stdout)
    except subprocess.CalledProcessError:
        print("samtools is not callable from everywhere. Please check your PATH.")
    except FileNotFoundError:
        print("samtools is not installed or not found in the PATH.")

def check_bcftools():
    try:
        # Check if bcftools is callable
        result = subprocess.run(['bcftools', '--version'], capture_output=True, text=True, check=True)
        print("checking bcftools")
        print(result.stdout)
    except subprocess.CalledProcessError:
        print("bcftools is not callable from everywhere. Please check your PATH.")
    except FileNotFoundError:
        print("bcftools is not installed or not found in the PATH.")

def check_agc():
    try:
        # Check if agc is callable
        result = subprocess.run(['agc', '--version'], capture_output=True, text=True, check=True)
        print("checking agc")
        print(result.stdout)
    except subprocess.CalledProcessError:
        print("agc is not callable from everywhere. Please check your PATH.")
    except FileNotFoundError:
        print("agc is not installed or not found in the PATH.")

def check_anchorwave():
    try:
        # Check if anchorwave is callable
        result = subprocess.run(['anchorwave', '--version'], capture_output=True, text=True, check=True)
        print("checking anchorwave")
        print(result.stdout)
    except subprocess.CalledProcessError:
        print("anchorwave is not callable from everywhere. Please check your PATH.")
    except FileNotFoundError:
        print("anchorwave is not installed or not found in the PATH.")

def check_python():
    try:
        # Check if python is callable
        result = subprocess.run(['python', '--version'], capture_output=True, text=True, check=True)
        print("checking python")
        print(result.stdout)
    except subprocess.CalledProcessError:
        print("python is not callable from everywhere. Please check your PATH.")
    except FileNotFoundError:
        print("python is not installed or not found in the PATH.")

def check_tiledb():   
    try:
        # Check if tiledb is callable
        result = subprocess.run(['tiledb', '--version'], capture_output=True, text=True, check=True)
        print("checking tiledb")
        print(result.stdout)
    except subprocess.CalledProcessError:
        print("tiledb is not callable from everywhere. Please check your PATH.")
    except FileNotFoundError:
        print("tiledb is not installed or not found in the PATH.")

def check_tiledbvcf():
    try:
        # Check if tiledbvcf is callable
        result = subprocess.run(['tiledbvcf', '--version'], capture_output=True, text=True, check=True)
        print("checking tiledbvcf")
        print(result.stdout)
    except subprocess.CalledProcessError:
        print("tiledbvcf is not callable from everywhere. Please check your PATH.")
    except FileNotFoundError:
        print("tiledbvcf is not installed or not found in the PATH.")

def check_perl():
    try:
        # Check if perl is callable
        result = subprocess.run(['perl', '--version'], capture_output=True, text=True, check=True)
        print("checking perl")
        print(result.stdout)
    except subprocess.CalledProcessError:
        print("perl is not callable from everywhere. Please check your PATH.")
    except FileNotFoundError:
        print("perl is not installed or not found in the PATH.")

def check_pygenometracks():
    try:
        # Check if pygenometracks is callable
        result = subprocess.run(['pyGenomeTracks', '--version'], capture_output=True, text=True, check=True)
        print("checking pygenometracks")
        print(result.stdout)
    except subprocess.CalledProcessError:
        print("pygenometracks is not callable from everywhere. Please check your PATH.")
    except FileNotFoundError:
        print("pygenometracks is not installed or not found in the PATH.")


def main():
    check_phgtools()
    check_phgv2()
    check_samtools()
    check_bcftools()
    check_agc()
    check_anchorwave()
    check_python()
    check_tiledb()
    check_tiledbvcf()
    check_perl()
    check_pygenometracks()


if __name__ == "__main__":
    main()