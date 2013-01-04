from pymol.wizard import Wizard
from pymol import cmd
from pika import *
import uuid
import threading
import json
import socket

WORDS_DICT = {'lys_linker': {'LYS': 'CA+CB+CB+CG+CD+CG'}, 'ser_linker': {'SER': 'CA+CB'}, 'alanine': {'ALA': 'CA+CB'}, 'guanidinium': {'ARG': 'CD+NE+CZ+NH1+CZ+NE+CZ+NH2'}, 'glu_linker': {'GLU': 'CA+CB+CB+CG+CD+CG'}, 'valine': {'VAL': 'CB+CG2+CA+CB+CB+CG1'}, 'lys_end': {'LYS': 'CD+CE+CE+NZ'}, 'isoleucine': {'ILE': 'CD1+CG1+CB+CG1+CB+CG2+CA+CB'}, 'thr_linker': {'THR': 'CA+CB+CB+CG2'}, 'indole': {'TRP': 'CD1+CG+CD1+NE1+CD2+CG+CE2+NE1+CD2+CE3+CD2+CE2+CE3+CZ3+CE2+CZ2+CH2+CZ3+CH2+CZ2'}, 'asn_linker': {'ASN': 'CB+CG+CA+CB'}, 'trp_linker': {'TRP': 'CA+CB+CB+CG'}, 'arg_linker': {'ARG': 'CD+CG+CB+CG+CA+CB'}, 'leucine': {'LEU': 'CA+CB+CB+CG+CD2+CG+CD1+CG'}, 'phospho_thr': {'TPO': 'O2P+P+OG+P+O1P+P+O3P+P+O1P+O3P'}, 'phospho_his': {'HIP': 'O2P+P+OG+P+O1P+P+O3P+P+O1P+O3P'}, 'proline': {'PRO': 'CD+CG+CB+CG'}, 'imidazole': {'HIS': 'CD2+NE2+CE1+NE2+CD2+CG+CE1+ND1+CG+ND1'}, 'carboxyl': {'ASP': 'CG+OD2+CG+OD1', 'GLU': 'CD+OE2+CD+OE1'}, 'phospho_tyr': {'PTR': 'O2P+P+OG+P+O1P+P+O3P+P+O1P+O3P'}, 'phenyl': {'PHE': 'CD2+CE2+CD2+CG+CE2+CZ+CD1+CG+CE1+CZ+CD1+CE1', 'TYR': 'CE1+CZ+CD1+CE1+CE2+CZ+CD1+CG+CD2+CE2+CD2+CG'}, 'asp_linker': {'ASP': 'CB+CG+CA+CB'}, 'hydroxyl': {'SER': 'CB+OG', 'TYR': 'CZ+OH', 'THR': 'CB+OG1'}, 'gln_linker': {'GLN': 'CA+CB+CB+CG+CD+CG'}, 'phospho_asp': {'PHD': 'O2P+P+OG+P+O1P+P+O3P+P+O1P+O3P'}, 'phospho_ser': {'SEP': 'O2P+P+OG+P+O1P+P+O3P+P+O1P+O3P'}, 'tyr_linker': {'TYR': 'CA+CB+CB+CG'}, 'peptide_bond': {'CYS': 'CA+N+C+CA+C+O', 'GLN': 'CA+N+C+CA+C+O', 'HIS': 'CA+N+C+CA+C+O', 'SER': 'CA+N+C+CA+C+O', 'VAL': 'CA+N+C+CA+C+O', 'LYS': 'CA+N+C+CA+C+O', 'ILE': 'CA+N+C+CA+C+O', 'PRO': 'CA+N+C+CA+C+O', 'GLY': 'CA+N+C+CA+C+O', 'THR': 'CA+N+C+CA+C+O', 'PHE': 'CA+N+C+CA+C+O', 'ALA': 'CA+N+C+CA+C+O', 'MET': 'CA+N+C+CA+C+O', 'ASP': 'CA+N+C+CA+C+O', 'GLU': 'CA+N+C+CA+C+O', 'LEU': 'CA+N+C+CA+C+O', 'ARG': 'CA+N+C+CA+C+O', 'TRP': 'CA+N+C+CA+C+O', 'ASN': 'CA+N+C+CA+C+O', 'TYR': 'CA+N+C+CA+C+O'}, 'carboxamide': {'ASN': 'CG+ND2+CG+OD1', 'GLN': 'CD+NE2+CD+OE1'}, 'phe_linker': {'PHE': 'CA+CB+CB+CG'}, 'cysteine': {'CYS': 'CA+CB+CB+SG'}, 'met_linker': {'MET': 'CA+CB+CB+CG'}, 'his_linker': {'HIS': 'CB+CG+CA+CB'}, 'met_end': {'MET': 'CE+SD+CG+SD'}}
BOND_WORD_DICT = {('PRO', 'CB', 'CG'): 'proline', ('ILE', 'C', 'O'): 'peptide_bond', ('TRP', 'CE2', 'NE1'): 'indole', ('ARG', 'CZ', 'NH1'): 'guanidinium', ('LYS', 'C', 'CA'): 'peptide_bond', ('TPO', 'O2P', 'P'): 'phospho_thr', ('HIP', 'O3P', 'P'): 'phospho_his', ('GLY', 'CA', 'N'): 'peptide_bond', ('ARG', 'CZ', 'NE'): 'guanidinium', ('ARG', 'C', 'CA'): 'peptide_bond', ('TYR', 'CA', 'N'): 'peptide_bond', ('TYR', 'CD2', 'CG'): 'phenyl', ('VAL', 'C', 'CA'): 'peptide_bond', ('SEP', 'O2P', 'P'): 'phospho_ser', ('TPO', 'O1P', 'O3P'): 'phospho_thr', ('GLY', 'C', 'CA'): 'peptide_bond', ('GLU', 'CD', 'CG'): 'glu_linker', ('TRP', 'CD1', 'NE1'): 'indole', ('TRP', 'CD2', 'CE2'): 'indole', ('ASP', 'C', 'CA'): 'peptide_bond', ('ILE', 'CB', 'CG1'): 'isoleucine', ('ASN', 'C', 'CA'): 'peptide_bond', ('MET', 'CB', 'CG'): 'met_linker', ('HIP', 'O2P', 'P'): 'phospho_his', ('PHE', 'CE1', 'CZ'): 'phenyl', ('PRO', 'C', 'CA'): 'peptide_bond', ('PHE', 'CE2', 'CZ'): 'phenyl', ('TYR', 'CA', 'CB'): 'tyr_linker', ('PHE', 'C', 'O'): 'peptide_bond', ('GLN', 'CD', 'NE2'): 'carboxamide', ('HIS', 'CG', 'ND1'): 'imidazole', ('SEP', 'O1P', 'P'): 'phospho_ser', ('PHE', 'CD1', 'CE1'): 'phenyl', ('PTR', 'O1P', 'O3P'): 'phospho_tyr', ('SEP', 'OG', 'P'): 'phospho_ser', ('ALA', 'CA', 'N'): 'peptide_bond', ('LEU', 'C', 'CA'): 'peptide_bond', ('THR', 'CB', 'OG1'): 'hydroxyl', ('ASP', 'CG', 'OD1'): 'carboxyl', ('PHD', 'O3P', 'P'): 'phospho_asp', ('LYS', 'CA', 'CB'): 'lys_linker', ('MET', 'C', 'O'): 'peptide_bond', ('TRP', 'CD2', 'CG'): 'indole', ('HIS', 'CB', 'CG'): 'his_linker', ('MET', 'CG', 'SD'): 'met_end', ('ARG', 'C', 'O'): 'peptide_bond', ('TRP', 'CA', 'CB'): 'trp_linker', ('ARG', 'CB', 'CG'): 'arg_linker', ('PTR', 'O3P', 'P'): 'phospho_tyr', ('LEU', 'CD1', 'CG'): 'leucine', ('SEP', 'O3P', 'P'): 'phospho_ser', ('MET', 'CA', 'CB'): 'met_linker', ('TYR', 'CD1', 'CE1'): 'phenyl', ('ALA', 'CA', 'CB'): 'alanine', ('GLU', 'CD', 'OE1'): 'carboxyl', ('PHE', 'CD1', 'CG'): 'phenyl', ('THR', 'CA', 'N'): 'peptide_bond', ('TYR', 'CZ', 'OH'): 'hydroxyl', ('ASP', 'CA', 'CB'): 'asp_linker', ('THR', 'CA', 'CB'): 'thr_linker', ('TRP', 'CA', 'N'): 'peptide_bond', ('SER', 'C', 'CA'): 'peptide_bond', ('ASP', 'CA', 'N'): 'peptide_bond', ('MET', 'CA', 'N'): 'peptide_bond', ('TRP', 'CH2', 'CZ2'): 'indole', ('ILE', 'CA', 'CB'): 'isoleucine', ('LYS', 'CD', 'CG'): 'lys_linker', ('PRO', 'CD', 'CG'): 'proline', ('GLU', 'C', 'O'): 'peptide_bond', ('TRP', 'CH2', 'CZ3'): 'indole', ('TPO', 'O1P', 'P'): 'phospho_thr', ('ILE', 'CA', 'N'): 'peptide_bond', ('CYS', 'C', 'CA'): 'peptide_bond', ('LYS', 'CD', 'CE'): 'lys_end', ('GLY', 'C', 'O'): 'peptide_bond', ('CYS', 'CA', 'CB'): 'cysteine', ('MET', 'C', 'CA'): 'peptide_bond', ('VAL', 'CA', 'N'): 'peptide_bond', ('THR', 'CB', 'CG2'): 'thr_linker', ('LEU', 'CA', 'CB'): 'leucine', ('ASP', 'C', 'O'): 'peptide_bond', ('LEU', 'CB', 'CG'): 'leucine', ('PRO', 'C', 'O'): 'peptide_bond', ('PHD', 'O1P', 'P'): 'phospho_asp', ('PHE', 'C', 'CA'): 'peptide_bond', ('ARG', 'CZ', 'NH2'): 'guanidinium', ('TYR', 'CB', 'CG'): 'tyr_linker', ('TYR', 'CE1', 'CZ'): 'phenyl', ('TYR', 'CE2', 'CZ'): 'phenyl', ('VAL', 'CB', 'CG2'): 'valine', ('LEU', 'C', 'O'): 'peptide_bond', ('LEU', 'CA', 'N'): 'peptide_bond', ('VAL', 'CA', 'CB'): 'valine', ('GLU', 'CA', 'N'): 'peptide_bond', ('LYS', 'C', 'O'): 'peptide_bond', ('SER', 'C', 'O'): 'peptide_bond', ('PTR', 'O2P', 'P'): 'phospho_tyr', ('GLU', 'CA', 'CB'): 'glu_linker', ('GLU', 'CB', 'CG'): 'glu_linker', ('GLN', 'CA', 'N'): 'peptide_bond', ('PHE', 'CD2', 'CG'): 'phenyl', ('TYR', 'C', 'CA'): 'peptide_bond', ('HIS', 'CA', 'CB'): 'his_linker', ('PTR', 'O1P', 'P'): 'phospho_tyr', ('TRP', 'CE2', 'CZ2'): 'indole', ('ASN', 'CG', 'OD1'): 'carboxamide', ('HIS', 'CA', 'N'): 'peptide_bond', ('ILE', 'CB', 'CG2'): 'isoleucine', ('HIS', 'CD2', 'NE2'): 'imidazole', ('TRP', 'CD2', 'CE3'): 'indole', ('ALA', 'C', 'CA'): 'peptide_bond', ('LEU', 'CD2', 'CG'): 'leucine', ('TRP', 'C', 'CA'): 'peptide_bond', ('GLN', 'CA', 'CB'): 'gln_linker', ('GLN', 'CB', 'CG'): 'gln_linker', ('HIP', 'OG', 'P'): 'phospho_his', ('THR', 'C', 'O'): 'peptide_bond', ('HIS', 'C', 'O'): 'peptide_bond', ('TRP', 'CD1', 'CG'): 'indole', ('ILE', 'C', 'CA'): 'peptide_bond', ('SER', 'CA', 'N'): 'peptide_bond', ('GLU', 'C', 'CA'): 'peptide_bond', ('GLN', 'C', 'CA'): 'peptide_bond', ('SEP', 'O1P', 'O3P'): 'phospho_ser', ('CYS', 'C', 'O'): 'peptide_bond', ('SER', 'CB', 'OG'): 'hydroxyl', ('HIS', 'CE1', 'NE2'): 'imidazole', ('PHE', 'CA', 'N'): 'peptide_bond', ('HIP', 'O1P', 'P'): 'phospho_his', ('ASP', 'CG', 'OD2'): 'carboxyl', ('ASN', 'CB', 'CG'): 'asn_linker', ('LYS', 'CE', 'NZ'): 'lys_end', ('CYS', 'CB', 'SG'): 'cysteine', ('VAL', 'CB', 'CG1'): 'valine', ('TYR', 'C', 'O'): 'peptide_bond', ('TPO', 'OG', 'P'): 'phospho_thr', ('VAL', 'C', 'O'): 'peptide_bond', ('HIP', 'O1P', 'O3P'): 'phospho_his', ('ARG', 'CD', 'NE'): 'guanidinium', ('ALA', 'C', 'O'): 'peptide_bond', ('LYS', 'CB', 'CG'): 'lys_linker', ('GLU', 'CD', 'OE2'): 'carboxyl', ('TPO', 'O3P', 'P'): 'phospho_thr', ('PHD', 'O2P', 'P'): 'phospho_asp', ('PHE', 'CA', 'CB'): 'phe_linker', ('ASN', 'CG', 'ND2'): 'carboxamide', ('HIS', 'CD2', 'CG'): 'imidazole', ('PHE', 'CD2', 'CE2'): 'phenyl', ('PRO', 'CA', 'N'): 'peptide_bond', ('ASN', 'C', 'O'): 'peptide_bond', ('ARG', 'CD', 'CG'): 'arg_linker', ('ARG', 'CA', 'N'): 'peptide_bond', ('TRP', 'CB', 'CG'): 'trp_linker', ('PHE', 'CB', 'CG'): 'phe_linker', ('HIS', 'C', 'CA'): 'peptide_bond', ('PHD', 'OG', 'P'): 'phospho_asp', ('MET', 'CE', 'SD'): 'met_end', ('GLN', 'C', 'O'): 'peptide_bond', ('HIS', 'CE1', 'ND1'): 'imidazole', ('CYS', 'CA', 'N'): 'peptide_bond', ('TRP', 'CE3', 'CZ3'): 'indole', ('TYR', 'CD2', 'CE2'): 'phenyl', ('ILE', 'CD1', 'CG1'): 'isoleucine', ('PTR', 'OG', 'P'): 'phospho_tyr', ('TRP', 'C', 'O'): 'peptide_bond', ('LYS', 'CA', 'N'): 'peptide_bond', ('TYR', 'CD1', 'CG'): 'phenyl', ('GLN', 'CD', 'CG'): 'gln_linker', ('ARG', 'CA', 'CB'): 'arg_linker', ('ASN', 'CA', 'N'): 'peptide_bond', ('SER', 'CA', 'CB'): 'ser_linker', ('GLN', 'CD', 'OE1'): 'carboxamide', ('PHD', 'O1P', 'O3P'): 'phospho_asp', ('THR', 'C', 'CA'): 'peptide_bond', ('ASN', 'CA', 'CB'): 'asn_linker', ('ASP', 'CB', 'CG'): 'asp_linker'}
SELECTION_NAME = 'query'
OBJECT_SUFFIX = 'suns'

SUNS_SERVER_ADDRESS = 'suns.degradolab.org'

################################################################################
# Here is the class that will actual do the search.
class SearchThread(threading.Thread):
    def __init__(self, rmsd, num_structures, random_seed, cmd, pdbstrs, server_address):
        threading.Thread.__init__(self)
        
        self.request = json.dumps({
            'rmsd_cutoff'   : rmsd,
            'num_structures': num_structures,
            'random_seed'   : random_seed,
            'atoms'         : pdbstrs
        })
        self.cmd = cmd
        self.channel = None
        self.callback_queue = None
        self.pdbs = {}
        self.current_status = "[*] Bug: 'current_status' unset"
        self.suns_server_address = server_address
        
    def handle_delivery(self, channel, method_frame, header_frame, body, corr_id=None):
        if(corr_id == header_frame.correlation_id):
            # TODO Consider changing to JSON in the future.
            # The first byte is the status.
            # 0 = No more search results
            # 1 = A search result
            # 2 = Search time limit exceeded
            # 3 = Error message
            if(body[0] == '0'):
                self.current_status = '[*] Search done.'
                self.channel.stop_consuming()
                return
            elif(body[0] == '2'):
                self.current_status = '[*] Time limit exceeded.'
                self.channel.stop_consuming()
                return
            elif(body[0] == '3'):
                print '[*] Error: ' + body[1:]
                return
                
            # The next 4 bytes are the pdbid.
            pdbid = body[1:5]
            if(pdbid not in self.pdbs):
                self.pdbs[pdbid] = 0
            sele_name = pdbid + '_%04d_%s' % (self.pdbs[pdbid], OBJECT_SUFFIX)
            # Delete this object if it already exists for some reason
            self.cmd.delete(sele_name)
            # Place the data in pymol.
            self.cmd.read_pdbstr(body[5:], sele_name)
            # Increment the number of results for this pdbid
            self.pdbs[pdbid] += 1
    
    # This will cancel any current searches.
    def cancel_search(self):
        # If we have declared a callback queue, then delete it
        # and stop consuming.
        if(self.callback_queue != None):
            self.current_status = '[*] Search cancelled.'
            self.channel.queue_delete(queue=self.callback_queue)
            self.channel.stop_consuming()
            
    # Iterate over the current results and delete them.
    def delete_current_results(self, exceptions = {}):
        for pdbid in self.pdbs:
            for i in range(self.pdbs[pdbid]):
                sele_name = pdbid + '_%04d_%s' % (i, OBJECT_SUFFIX)
                if(sele_name not in exceptions):
                    self.cmd.delete(sele_name)
    
    def get_current_results(self):
        return self.pdbs
    
    # Perform our search.
    def run(self):
        # Create a unique identifier.
        try:
            corr_id = str(uuid.uuid4())
            ############################################################################
            # First, setup the channel to the server.
            credentials = PlainCredentials('suns-client', 'suns-client')
            connection = BlockingConnection(ConnectionParameters(host=self.suns_server_address, credentials=credentials, virtual_host='suns-vhost'))
            self.channel = connection.channel()
            self.channel.exchange_declare(exchange='suns-exchange-responses', passive=True, durable=True)
            self.channel.exchange_declare(exchange='suns-exchange-requests', passive=True, durable=True)
            # Now create the callback queue
            result = self.channel.queue_declare(exclusive=True)
            self.callback_queue = result.method.queue
            self.channel.queue_bind(exchange='suns-exchange-responses', queue = self.callback_queue, routing_key=self.callback_queue)
            self.channel.basic_consume(lambda c, m, h, b : self.handle_delivery(c, m, h, b, corr_id), no_ack=True, queue=self.callback_queue)
            ############################################################################
            # Now ask the server to perform the search.
            print '[*] Searching...'
            # Now issue the request.
            self.channel.basic_publish(exchange = 'suns-exchange-requests', routing_key = '1.0.0',
                                  properties=BasicProperties(reply_to = self.callback_queue, correlation_id = corr_id ),
                                  body = self.request)

            # This call is a blocking call.  It will block until channel.stop_consuming is called.
            self.channel.start_consuming()

            self.cmd.orient(SELECTION_NAME)
            print self.current_status

            ############################################################################
            # Now close the connection
            connection.close()
        except socket.error:
            print "[*] Error: Unable to connect to '" + self.suns_server_address + "'"

################################################################################
# Wizard class
class Suns_search(Wizard):
    '''
    This class  will create the wizard for performing Suns searches.
    '''
    def __init__(self, _self=cmd):
        Wizard.__init__(self, _self)
        self.cmd.unpick()
        self.word_list = {}
        self.current_selection = ''
        self.rmsd_cutoff = 1.0
        self.random_seed = 0
        self.number_of_structures = 100
        self.searchThread = None
        self.prev_auto_hide_setting = self.cmd.get('auto_hide_selection')
        cmd.set('auto_hide_selections',0)
        self.prev_mouse_mode = self.cmd.get('mouse_selection_mode')
        cmd.set('mouse_selection_mode',0)
        cmd.config_mouse('three_button_editing')
        self.do_select(SELECTION_NAME, 'none')
        self.suns_server_address = SUNS_SERVER_ADDRESS
    
    def cleanup(self):
        '''
        Once we are done with the wizard, we should set various pymol
        parameters back to their original values.
        '''
        cmd.config_mouse('three_button_viewing')
        cmd.set('auto_hide_selections', self.prev_auto_hide_setting)
        cmd.set('mouse_selection_mode', self.prev_mouse_mode)
        self.cmd.delete(SELECTION_NAME)
    
    def set_rmsd(self, rmsd):
        '''
        This is the method that will be called once the user has
        selected an rmsd cutoff via the wizard menu.
        '''
        self.rmsd_cutoff = rmsd
        self.cmd.refresh_wizard()
    
    def set_num_structures(self, num_structures):
        '''
        This is the method that will be called once the user
        has set the maximum number of structures to return.
        '''
        self.number_of_structures = num_structures
        self.cmd.refresh_wizard()
    
    def get_random_seed(self):
        '''
        This method will bring up a dialog box and prompt the user for
        an integer to use as the random seed in shuffling the Suns database.
        '''
        import Tkinter
        import tkSimpleDialog
        Tkinter.Tk()
        random_seed = tkSimpleDialog.askstring('Randomize Order','Random Seed (Integer):')
        return random_seed
    
    def ask_random_seed(self, ask_for_random_seed):
        ''' 
        This method will be called when the user clicks on the wizard
        menu to randomize the Suns database.  If the user selects
        to randomize the db, then this method will call get_random_seed
        to get the random seed for randomization.
        '''
        if(ask_for_random_seed):
            try:
                self.random_seed = int(self.get_random_seed())
            except:
                self.random_seed = 0
        else:
            self.random_seed = 0
        self.cmd.refresh_wizard()
    
    def create_random_seed_menu(self):
        '''
        This method will create the Randomize Database menu.
        '''
        random_seed_menu = [[2, 'Randomize Order', '']]
        random_seed_menu.append([1, 'No', 'cmd.get_wizard().ask_random_seed(False)'])
        random_seed_menu.append([1, 'Yes', 'cmd.get_wizard().ask_random_seed(True)'])
        
        return random_seed_menu
        
    def create_rmsd_menu(self):
        '''
        This method will create a wizard menu for the possible RMSD cutoff values.
        Currently the values range from 0.1 to 2 A RMSD.
        '''
        rmsd_menu = [[2, 'RMSD Cutoff', '']]
        for rmsd_choice in range(1,21):
            rmsd = float(rmsd_choice) / 10.0
            rmsd_menu.append([1, str(rmsd) , 'cmd.get_wizard().set_rmsd(' + str(rmsd) + ')'])
        return rmsd_menu

    def create_num_structures_menu(self):
        '''
        This method will create a wizard menu for the possible number of structures
        to return.  Values range from 10 to 2000.
        '''
        num_structures_menu = [[2, 'Number of Results', '']]
        num_structures_menu.append([1, str(10) , 'cmd.get_wizard().set_num_structures(' + str(10) + ')'])
        num_structures_menu.append([1, str(50) , 'cmd.get_wizard().set_num_structures(' + str(50) + ')'])
        num_structures_menu.append([1, str(100) , 'cmd.get_wizard().set_num_structures(' + str(100) + ')'])
        num_structures_menu.append([1, str(250) , 'cmd.get_wizard().set_num_structures(' + str(250) + ')'])
        num_structures_menu.append([1, str(500) , 'cmd.get_wizard().set_num_structures(' + str(500) + ')'])
        num_structures_menu.append([1, str(1000) , 'cmd.get_wizard().set_num_structures(' + str(1000) + ')'])
        num_structures_menu.append([1, str(1500) , 'cmd.get_wizard().set_num_structures(' + str(1500) + ')'])
        num_structures_menu.append([1, str(2000) , 'cmd.get_wizard().set_num_structures(' + str(2000) + ')'])
        
        return num_structures_menu
    
    def create_server_address_menu(self):
        '''
        This method will create the Server Address menu.
        '''
        server_address_menu = [[2, 'Server Address', '']]
        server_address_menu.append([1, 'Default: ' + SUNS_SERVER_ADDRESS, 'cmd.get_wizard().set_server_address("' + SUNS_SERVER_ADDRESS + '")'])
        server_address_menu.append([1, 'User defined', 'cmd.get_wizard().ask_server_address()'])
        
        return server_address_menu
        
    def set_server_address(self, server_address):
        '''
        This method will set the server address to use.
        '''
        self.suns_server_address = server_address
        self.cmd.refresh_wizard()
        
    def ask_server_address(self):
        '''
        This method will bring up a dialog box and prompt the user for
        the address to the Suns server.  This will allow the user to specify
        an address.
        '''
        import Tkinter
        import tkSimpleDialog
        Tkinter.Tk()
        server_address = tkSimpleDialog.askstring('Suns server address','Suns server address:')
        
        if server_address: self.set_server_address(server_address)
    
    def get_panel(self):
        '''
        This is the wizard method that will create the main menu.
        '''
        rmsd_menu = self.create_rmsd_menu()
        self.menu['rmsd'] = rmsd_menu
        num_structures_menu = self.create_num_structures_menu()
        self.menu['num_structures'] = num_structures_menu
        random_seed_menu = self.create_random_seed_menu()
        self.menu['random_seed'] = random_seed_menu
        server_address_menu = self.create_server_address_menu()
        self.menu['server'] = server_address_menu
        
        return [
            [ 1, 'Structural Search Engine',''],
            [ 2, 'Search', 'cmd.get_wizard().launch_search()'],
            [ 2, 'Cancel Search','cmd.get_wizard().cancel_search()'],
            [ 3, 'RMSD Cutoff: ' + str(self.rmsd_cutoff) + ' Angstroms', 'rmsd'],
            [ 3, 'Cap: ' + str(self.number_of_structures) + ' results', 'num_structures'],
            [ 3, 'Order: ' + {True: 'Random (Seed = %d' % self.random_seed + ')', False: 'Default'}[self.random_seed != 0], 'random_seed'],
            [ 3, 'Server: ' + self.suns_server_address, 'server'],
            [ 2, 'Clear Results', 'cmd.get_wizard().delete_current_results()'],
            [ 2, 'Clear Selection','cmd.get_wizard().clear_selection()'],
            [ 2, 'Fetch Full Context','cmd.get_wizard().fetch_full_context()'],
# There should be a cancel for Fetch Full Context
            [ 2, 'Done','cmd.set_wizard()'] ]
    
    def fetch_full_context(self):
        '''
        This method will look at what search results are currently enabled,
        and then fetch and align the full structure so that the user can 
        see the full context of the match.
        '''
        if(self.searchThread == None):
            return
        
        #current_results = self.searchThread.get_current_results()
        objs = cmd.get_names('objects', enabled_only=1)
        for obj in objs:
            w = obj.split('_')
            # Note that the search results will be named pdbid_XXXX.
            if( (len(w) == 2) and (len(w[0]) == 4) and (len(w[1]) == 4)):
                new_object_name = w[0] + '_full_' + w[1]
                cmd.fetch(w[0], new_object_name)
                dict = {'x' : []}
                # Get the atom info for the current object.
                self.cmd.iterate(obj, 'x.append( [model,segi,chain,resn,resi,name,alt] )', space=dict)
                selection = []
                for item in dict['x']:
                    selection += ['(chain %s and resn %s and resi %s and name %s)' % tuple(item[2:6])]
                selection = '(' + ' or '.join(selection) + ')'
                cmd.super(new_object_name + ' and (' + selection + ')', obj + ' and (' + selection + ')')
                
    
    def delete_current_results(self, exceptions={}):
        '''
        This method will ask the search thread (who knows all of the current
        results) to delete all results with the exception of any objects
        that are in the exceptions variable.
        '''
        if(self.searchThread != None):
            self.searchThread.delete_current_results(exceptions)
    
    def get_current_object_names(self, selectionName):
        '''
        This method will get all of the objects in the given selection.
        It is used so that when we perform a new search, we first check to see
        what objects are part of the search, and don't delete them.
        '''
        currentObjects = {}
        dict = {'x' : []}
        # Get the atom info for the two atoms of the bond.
        self.cmd.iterate(selectionName, 'x.append( (model,segi,chain,resn,resi,name,alt) )', space=dict)
        for a in dict['x']:
            currentObjects[a[0]] = 1
        
        return currentObjects
    
    def launch_search(self):
        '''
        This method will actually launch the search.
        '''
        pdbstr = self.cmd.get_pdbstr(SELECTION_NAME)
        exceptions = self.get_current_object_names(SELECTION_NAME)
        if(self.searchThread != None):
            self.searchThread.cancel_search()
            self.delete_current_results(exceptions)
        self.searchThread = SearchThread(self.rmsd_cutoff, self.number_of_structures, self.random_seed, self.cmd, pdbstr, self.suns_server_address)
        self.searchThread.start()
    
    def clear_selection(self):
        '''
        This method simply deletes the current selection.
        '''
        self.word_list = {}
        self.cmd.delete(SELECTION_NAME)
    
    def cancel_search(self):
        '''
        This method will cancel searching.
        '''
        if(self.searchThread != None):
            self.searchThread.cancel_search()
    
    def do_select(self, name, selection_logic):
        if(selection_logic != ''):
            self.cmd.select(name, selection_logic)
            self.cmd.enable(name)
            
    def do_pick(self, bondFlag):
        '''
        This is the method that is called each time the user
        uses the mouse to select a bond or atom.
        '''
        # I think bondFlag only = 1 if we are selecting bonds.
        # We are only accepting bond selections.
        if(bondFlag == 1):
            dict = {'x' : []}
            # Get the atom info for the two atoms of the bond.
            self.cmd.iterate("pk1", 'x.append( (model,segi,chain,resn,resi,name,alt) )', space=dict)
            self.cmd.iterate("pk2", 'x.append( (model,segi,chain,resn,resi,name,alt) )', space=dict)
            bond = [dict['x'][0][5], dict['x'][1][5]]
            bond.sort()
            
            # We only handle bonds per residue, no residue spanning bonds.
            if( (dict['x'][0][2] != dict['x'][1][2]) or (dict['x'][0][4] != dict['x'][1][4]) ):
                return
            
            # The key is of the form (resn, atom0, atom1) where the atom names are in alphabetical order.
            key = (dict['x'][0][3], bond[0], bond[1])
            # Look up what word this bond is part of.
            if(key in BOND_WORD_DICT):
                word = BOND_WORD_DICT[key]
                # Now form a new key that has the current residue and the word.
                key = tuple(list(dict['x'][0][0:5]) + [word])
                if(key in self.word_list):
                    del self.word_list[key]
                else:
                    if(key[1].strip() == ''):
                        self.word_list[key] = '(model %s and chain %s and resn %s and resi %s and name %s )' % tuple(list(key[0:1]) + list(key[2:5]) + [WORDS_DICT[word][key[3]]])
                    else:
                        self.word_list[key] = '(model %s and segi %s and chain %s and resn %s and resi %s and name %s )' % tuple(list(key[0:5]) + [WORDS_DICT[word][key[3]]])
            #else: # If this isn't part of a word, we'll still let them select the bond.
            #    # If this bond isn't part of a word, then the key
            #    # will include the full identifier for the two atoms.
            #    # We should sort the atoms so that if the bond comes in as a,b or b,a
            #    # we'll know it's the same thing.  For now, I'll just sort on the atom name.
            #    # This could cause a problem if the two atoms have the same name, but really,
            #    # should a bond have two atoms of the same name?
            #    if(dict['x'][0][5] < dict['x'][1][5]):
            #        key = tuple(list(dict['x'][0][0:6]) + list(dict['x'][1][0:6]))
            #    else:
            #        key = tuple(list(dict['x'][1][0:6]) + list(dict['x'][0][0:6]))
            #    if(key in self.word_list):
            #        del self.word_list[key]
            #    else:
            #        self.word_list[key] = '((model %s and segi %s and chain %s and resn %s and resi %s and name %s) or (model %s and segi %s and chain %s and resn %s and resi %s and name %s) )' % key
            
            self.cmd.unpick()
            self.current_selection = ' or '.join(self.word_list.values()).strip()
            self.do_select(SELECTION_NAME, self.current_selection)
