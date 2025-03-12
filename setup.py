import os
import subprocess
import sys

def setup_backend():
    """Set up the backend environment and install dependencies."""
    print("Setting up backend...")
    os.chdir('backend')
    
    # Create virtual environment
    if not os.path.exists('venv'):
        print("Creating virtual environment...")
        if sys.platform == 'win32':
            subprocess.run([sys.executable, '-m', 'venv', 'venv'])
        else:
            subprocess.run([sys.executable, '-m', 'venv', 'venv'])
    
    # Activate virtual environment and install dependencies
    print("Installing backend dependencies...")
    if sys.platform == 'win32':
        pip_path = os.path.join('venv', 'Scripts', 'pip')
        subprocess.run([pip_path, 'install', '-r', 'requirements.txt'])
    else:
        pip_path = os.path.join('venv', 'bin', 'pip')
        subprocess.run([pip_path, 'install', '-r', 'requirements.txt'])
    
    os.chdir('..')

def setup_frontend():
    """Set up the frontend and install dependencies."""
    print("Setting up frontend...")
    os.chdir('frontend')
    
    # Install npm dependencies
    print("Installing frontend dependencies...")
    if sys.platform == 'win32':
        subprocess.run(['npm', 'install'])
    else:
        subprocess.run(['npm', 'install'])
    
    os.chdir('..')

def main():
    """Set up both the frontend and backend."""
    print("Setting up Drug Discovery Assistant...")
    
    # Set up backend
    setup_backend()
    
    # Set up frontend
    setup_frontend()
    
    print("\nSetup complete!")
    print("\nTo run the application, use:")
    print("python run.py")

if __name__ == "__main__":
    main() 