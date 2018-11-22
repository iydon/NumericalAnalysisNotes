import platform
import os


system = platform.system()
ignore = ["aux", "log", "synctex", "thm", "toc", "out"]

if system=="Windows":
	cmd = ["del *.%s"%i for i in ignore]
	os.system(" ; ".join(cmd))
elif system=="Linux":
	cmd = ["rm *.%s"%i for i in ignore]
	os.system(" ; ".join(cmd))
else:
	print("Not Support.")
