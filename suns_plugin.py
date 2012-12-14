from pymol import cmd

def __init__(self):
    self.menuBar.addmenuitem(
        'Plugin',
        'command',
        label = 'Suns Search',
        command = lambda s=self:search() )

    def search():
        cmd.wizard('suns_search')
