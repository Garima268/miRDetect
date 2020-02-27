import importlib.util
import os
import sys
import subprocess

package_name = 'pyfasta'
ML_name = 'pandas'

#Check whether Pyfasta is installed or not
spec = importlib.util.find_spec(package_name)
if spec is None:
    print(package_name +" is not installed")
    print("Installing" + package_name)
    subprocess.run(["pip", "install", "pyfasta"])
else:
    print(package_name + " is installed")
  
#Check whether Sklearn is installed or not
spec = importlib.util.find_spec(ML_name)
if spec is None:
    print(ML_name +" is not installed")
    print("Installing" + ML_name)
    subprocess.run(["pip", "install", "pandas"])
else:
    print(ML_name + " is installed")
  