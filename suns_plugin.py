from pymol import cmd

def __init__(self):
    self.menuBar.addmenuitem(
        'Wizard',
        'command',
        label = 'Search',
        command = lambda s=self:search() )

    def search():
        cmd.wizard('suns_search')
