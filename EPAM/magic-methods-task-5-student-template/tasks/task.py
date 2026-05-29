import os
import shutil
from uuid import uuid4


class TempDir:
    def __enter__(self):
        self.current_dir = os.getcwd()
        print(f"Current dir path: {os.getcwd()}")
        dirname = str(uuid4())

        os.mkdir(os.path.join(os.getcwd(),dirname))
        self.path = os.path.join(os.getcwd(),dirname)
        os.chdir(self.path)
        print(f"New directory path at: {os.getcwd()}")
        return f"Current new path: {self.path}"
    
    def __exit__(self,exc_type,exc_val,exc_tb):
        print(f"Removing current dir: {self.path}")
        os.chdir(self.current_dir)
        shutil.rmtree(self.path)
        print(f"Current dir is now: {os.getcwd()}")
       

with TempDir() as result:
    print(f"This is result string: {result}")