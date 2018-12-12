
import pathlib
import os

class FileOpen:

    def __init__(self, folder, filename):
        folder_path = pathlib.Path(folder)
        if not folder_path.is_dir():
            raise NotADirectoryError("Folder " + str(folder) + " does not exist.")
            # 'folder' must be enclosed within str() otherwise: TypeError: must be str, not WindowsPath

        self.absolute = pathlib.Path(os.path.join(folder, filename))
        if not self.absolute.exists():
            raise FileNotFoundError("File +" + str(self.absolute) + " does not exist.")

        # https://docs.python.org/3/library/os.path.html#os.path.splitext
        self.filename, self.postfix = os.path.splitext(filename)
        self.folder = folder

    def __str__(self):
        return "FileOpen:\n" \
               + "Filename: " + str(self.filename) + "\n" \
               + "Postfix: " + str(self.postfix) + "\n" \
               + "Folder: " + str(self.folder) + "\n" \
               + "Absolute: " + str(self.absolute)
