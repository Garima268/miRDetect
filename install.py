import importlib.util
import subprocess
import sys


def ensure_package(package_name, install_name=None):
    """
    Check if a package is installed; if not, install it.
    
    Args:
        package_name (str): name used for import
        install_name (str): name used for pip install (if different)
    """
    install_name = install_name or package_name

    if importlib.util.find_spec(package_name) is None:
        print(f"{package_name} is not installed. Installing...")
        try:
            subprocess.check_call([sys.executable, "-m", "pip", "install", install_name])
            print(f"{install_name} installed successfully.")
        except subprocess.CalledProcessError:
            print(f"Failed to install {install_name}.")
    else:
        print(f"{package_name} is already installed.")


def main():
    packages = [
        ("pyfasta", None),
        ("pandas", None),
    ]

    for pkg, install_name in packages:
        ensure_package(pkg, install_name)


if __name__ == "__main__":
    main()
