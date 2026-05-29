import os
import uuid
from shutil import rmtree


class TempDir:
    def __enter__(self):
        self.previous_dir = os.getcwd()

        self.temp_dir = os.path.join(
            self.previous_dir,
            str(uuid.uuid4())
        )

        os.mkdir(self.temp_dir)
        os.chdir(self.temp_dir)

        return self.temp_dir

    def __exit__(self, exc_type, exc_val, exc_tb):
        os.chdir(self.previous_dir)
        rmtree(self.temp_dir)