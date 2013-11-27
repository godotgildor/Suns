import sys
import base64

INSTALL_CODE = '''import os
import base64
import zipfile
import tempfile

class SunsInstaller:
    def __init__(self):
        dir = os.path.dirname(os.path.abspath(__file__))
        tmpZip = tempfile.NamedTemporaryFile(delete=False)
        tmpZip.write(base64.b64decode(CODE))
        tmpZip.close()

        with zipfile.ZipFile(tmpZip.name, 'r') as zf:
            zf.extractall(dir)

        os.unlink(tmpZip.name)

def __init__(self):
    self.menuBar.addmenuitem('Plugin', 'command', 'Suns Installer',
                             label = 'Suns Installer',
                             command = lambda s = self : SunsInstaller(s))
'''


if __name__ == '__main__':
    if(len(sys.argv) < 2):
        print 'Usage: python make_py_installer.py <zip installer> <py installer>'
        sys.exit(1)

    txt = base64.b64encode(open(sys.argv[1], 'rb').read())
    of = open(sys.argv[2], 'w')
    of.write('CODE = """' + txt + '"""\n\n')
    of.write(INSTALL_CODE)
    of.close()
