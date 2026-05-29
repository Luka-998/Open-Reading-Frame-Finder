import os
import shutil

class TempDir:
    def __enter__(self):
        current_dir = os.getcwd()
        self.current_dir = current_dir
        print(f"Current dir path: {os.getcwd()}")
        os.mkdir(os.path.join(os.getcwd(),'new_dir'))
        path = os.path.join(os.getcwd(),'new_dir')
        os.chdir(path)
        print(f"New directory path at: {os.getcwd()}")
       
    def __exit__(self,exc_type,exc_val,exc_tb):
        print(f"Removing current dir: {os.getcwd()}")
        os.chdir(self.current_dir)
        shutil.rmtree(os.path.join(os.getcwd(),'new_dir'))
        print(f"Current dir is now: {os.getcwd()}")
       

with TempDir():
    print("Yes")