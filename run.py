import os
import subprocess
import sys
import time
from threading import Thread

def run_backend():
    """Run the FastAPI backend server."""
    os.chdir('backend')
    if sys.platform == 'win32':
        venv_python = os.path.join('venv', 'Scripts', 'python')
        subprocess.run([venv_python, '-m', 'uvicorn', 'main:app', '--reload', '--host', '0.0.0.0', '--port', '8000'])
    else:
        venv_python = os.path.join('venv', 'bin', 'python')
        subprocess.run([venv_python, '-m', 'uvicorn', 'main:app', '--reload', '--host', '0.0.0.0', '--port', '8000'])

def run_frontend():
    """Run the React frontend development server."""
    os.chdir('frontend')
    if sys.platform == 'win32':
        subprocess.run(['npm', 'start'])
    else:
        subprocess.run(['npm', 'start'])

def main():
    """Run both the frontend and backend servers."""
    print("Starting Drug Discovery Assistant...")
    
    # Start the backend in a separate thread
    backend_thread = Thread(target=run_backend)
    backend_thread.daemon = True
    backend_thread.start()
    
    # Give the backend a moment to start
    print("Starting backend server...")
    time.sleep(2)
    
    # Start the frontend
    print("Starting frontend server...")
    run_frontend()

if __name__ == "__main__":
    main() 