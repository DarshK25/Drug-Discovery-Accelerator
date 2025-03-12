import os
import subprocess
import sys

def main():
    """Install backend dependencies."""
    print("Installing backend dependencies...")
    
    # Change to backend directory
    os.chdir('backend')
    
    # Create virtual environment if it doesn't exist
    if not os.path.exists('venv'):
        print("Creating virtual environment...")
        subprocess.run([sys.executable, '-m', 'venv', 'venv'])
    
    # Install dependencies
    print("Installing packages...")
    if sys.platform == 'win32':
        pip_path = os.path.join('venv', 'Scripts', 'pip')
    else:
        pip_path = os.path.join('venv', 'bin', 'pip')
    
    subprocess.run([pip_path, 'install', '-r', 'requirements.txt'])
    
    print("Backend dependencies installed successfully!")

if __name__ == "__main__":
    main() 