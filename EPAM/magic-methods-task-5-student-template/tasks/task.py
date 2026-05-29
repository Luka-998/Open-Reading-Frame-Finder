import os
import shutil
import uuid


class TempDir:
    def __enter__(self):
        self.previous_dir = os.getcwd()

        dirname = str(uuid.uuid4())
        self.temp_dir = os.path.join(self.previous_dir, dirname)

        os.mkdir(self.temp_dir)
        os.chdir(self.temp_dir)

        return self.temp_dir

    def __exit__(self, exc_type, exc_val, exc_tb):
        os.chdir(self.previous_dir)
        shutil.rmtree(self.temp_dir)