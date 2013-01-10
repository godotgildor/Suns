from pymol.wizard import Wizard
from pymol import cmd
from pika import *
import pymol.controlling
import uuid
import threading
import json
import socket
import Tkinter
import tkSimpleDialog
import tkMessageBox

# The following two dictionaries are automatically generated from the motif
# directories
WORDS_DICT = {'lys_linker': {'LYS': 'CA+CB+CB+CG+CD+CG'}, 'ser_linker': {'SER': 'CA+CB'}, 'alanine': {'ALA': 'CA+CB'}, 'guanidinium': {'ARG': 'CD+NE+CZ+NH1+CZ+NE+CZ+NH2'}, 'glu_linker': {'GLU': 'CA+CB+CB+CG+CD+CG'}, 'valine': {'VAL': 'CB+CG2+CA+CB+CB+CG1'}, 'lys_end': {'LYS': 'CD+CE+CE+NZ'}, 'isoleucine': {'ILE': 'CD1+CG1+CB+CG1+CB+CG2+CA+CB'}, 'thr_linker': {'THR': 'CA+CB+CB+CG2'}, 'indole': {'TRP': 'CD1+CG+CD1+NE1+CD2+CG+CE2+NE1+CD2+CE3+CD2+CE2+CE3+CZ3+CE2+CZ2+CH2+CZ3+CH2+CZ2'}, 'asn_linker': {'ASN': 'CB+CG+CA+CB'}, 'trp_linker': {'TRP': 'CA+CB+CB+CG'}, 'arg_linker': {'ARG': 'CD+CG+CB+CG+CA+CB'}, 'leucine': {'LEU': 'CA+CB+CB+CG+CD2+CG+CD1+CG'}, 'phospho_thr': {'TPO': 'O2P+P+OG+P+O1P+P+O3P+P+O1P+O3P'}, 'phospho_his': {'HIP': 'O2P+P+OG+P+O1P+P+O3P+P+O1P+O3P'}, 'proline': {'PRO': 'CD+CG+CB+CG'}, 'imidazole': {'HIS': 'CD2+NE2+CE1+NE2+CD2+CG+CE1+ND1+CG+ND1'}, 'carboxyl': {'ASP': 'CG+OD2+CG+OD1', 'GLU': 'CD+OE2+CD+OE1'}, 'phospho_tyr': {'PTR': 'O2P+P+OG+P+O1P+P+O3P+P+O1P+O3P'}, 'phenyl': {'PHE': 'CD2+CE2+CD2+CG+CE2+CZ+CD1+CG+CE1+CZ+CD1+CE1', 'TYR': 'CE1+CZ+CD1+CE1+CE2+CZ+CD1+CG+CD2+CE2+CD2+CG'}, 'asp_linker': {'ASP': 'CB+CG+CA+CB'}, 'hydroxyl': {'SER': 'CB+OG', 'TYR': 'CZ+OH', 'THR': 'CB+OG1'}, 'gln_linker': {'GLN': 'CA+CB+CB+CG+CD+CG'}, 'phospho_asp': {'PHD': 'O2P+P+OG+P+O1P+P+O3P+P+O1P+O3P'}, 'phospho_ser': {'SEP': 'O2P+P+OG+P+O1P+P+O3P+P+O1P+O3P'}, 'tyr_linker': {'TYR': 'CA+CB+CB+CG'}, 'peptide_bond': {'CYS': 'CA+N+C+CA+C+O', 'GLN': 'CA+N+C+CA+C+O', 'HIS': 'CA+N+C+CA+C+O', 'SER': 'CA+N+C+CA+C+O', 'VAL': 'CA+N+C+CA+C+O', 'LYS': 'CA+N+C+CA+C+O', 'ILE': 'CA+N+C+CA+C+O', 'PRO': 'CA+N+C+CA+C+O', 'GLY': 'CA+N+C+CA+C+O', 'THR': 'CA+N+C+CA+C+O', 'PHE': 'CA+N+C+CA+C+O', 'ALA': 'CA+N+C+CA+C+O', 'MET': 'CA+N+C+CA+C+O', 'ASP': 'CA+N+C+CA+C+O', 'GLU': 'CA+N+C+CA+C+O', 'LEU': 'CA+N+C+CA+C+O', 'ARG': 'CA+N+C+CA+C+O', 'TRP': 'CA+N+C+CA+C+O', 'ASN': 'CA+N+C+CA+C+O', 'TYR': 'CA+N+C+CA+C+O'}, 'carboxamide': {'ASN': 'CG+ND2+CG+OD1', 'GLN': 'CD+NE2+CD+OE1'}, 'phe_linker': {'PHE': 'CA+CB+CB+CG'}, 'cysteine': {'CYS': 'CA+CB+CB+SG'}, 'met_linker': {'MET': 'CA+CB+CB+CG'}, 'his_linker': {'HIS': 'CB+CG+CA+CB'}, 'met_end': {'MET': 'CE+SD+CG+SD'}}
BOND_WORD_DICT = {('PRO', 'CB', 'CG'): 'proline', ('ILE', 'C', 'O'): 'peptide_bond', ('TRP', 'CE2', 'NE1'): 'indole', ('ARG', 'CZ', 'NH1'): 'guanidinium', ('LYS', 'C', 'CA'): 'peptide_bond', ('TPO', 'O2P', 'P'): 'phospho_thr', ('HIP', 'O3P', 'P'): 'phospho_his', ('GLY', 'CA', 'N'): 'peptide_bond', ('ARG', 'CZ', 'NE'): 'guanidinium', ('ARG', 'C', 'CA'): 'peptide_bond', ('TYR', 'CA', 'N'): 'peptide_bond', ('TYR', 'CD2', 'CG'): 'phenyl', ('VAL', 'C', 'CA'): 'peptide_bond', ('SEP', 'O2P', 'P'): 'phospho_ser', ('TPO', 'O1P', 'O3P'): 'phospho_thr', ('GLY', 'C', 'CA'): 'peptide_bond', ('GLU', 'CD', 'CG'): 'glu_linker', ('TRP', 'CD1', 'NE1'): 'indole', ('TRP', 'CD2', 'CE2'): 'indole', ('ASP', 'C', 'CA'): 'peptide_bond', ('ILE', 'CB', 'CG1'): 'isoleucine', ('ASN', 'C', 'CA'): 'peptide_bond', ('MET', 'CB', 'CG'): 'met_linker', ('HIP', 'O2P', 'P'): 'phospho_his', ('PHE', 'CE1', 'CZ'): 'phenyl', ('PRO', 'C', 'CA'): 'peptide_bond', ('PHE', 'CE2', 'CZ'): 'phenyl', ('TYR', 'CA', 'CB'): 'tyr_linker', ('PHE', 'C', 'O'): 'peptide_bond', ('GLN', 'CD', 'NE2'): 'carboxamide', ('HIS', 'CG', 'ND1'): 'imidazole', ('SEP', 'O1P', 'P'): 'phospho_ser', ('PHE', 'CD1', 'CE1'): 'phenyl', ('PTR', 'O1P', 'O3P'): 'phospho_tyr', ('SEP', 'OG', 'P'): 'phospho_ser', ('ALA', 'CA', 'N'): 'peptide_bond', ('LEU', 'C', 'CA'): 'peptide_bond', ('THR', 'CB', 'OG1'): 'hydroxyl', ('ASP', 'CG', 'OD1'): 'carboxyl', ('PHD', 'O3P', 'P'): 'phospho_asp', ('LYS', 'CA', 'CB'): 'lys_linker', ('MET', 'C', 'O'): 'peptide_bond', ('TRP', 'CD2', 'CG'): 'indole', ('HIS', 'CB', 'CG'): 'his_linker', ('MET', 'CG', 'SD'): 'met_end', ('ARG', 'C', 'O'): 'peptide_bond', ('TRP', 'CA', 'CB'): 'trp_linker', ('ARG', 'CB', 'CG'): 'arg_linker', ('PTR', 'O3P', 'P'): 'phospho_tyr', ('LEU', 'CD1', 'CG'): 'leucine', ('SEP', 'O3P', 'P'): 'phospho_ser', ('MET', 'CA', 'CB'): 'met_linker', ('TYR', 'CD1', 'CE1'): 'phenyl', ('ALA', 'CA', 'CB'): 'alanine', ('GLU', 'CD', 'OE1'): 'carboxyl', ('PHE', 'CD1', 'CG'): 'phenyl', ('THR', 'CA', 'N'): 'peptide_bond', ('TYR', 'CZ', 'OH'): 'hydroxyl', ('ASP', 'CA', 'CB'): 'asp_linker', ('THR', 'CA', 'CB'): 'thr_linker', ('TRP', 'CA', 'N'): 'peptide_bond', ('SER', 'C', 'CA'): 'peptide_bond', ('ASP', 'CA', 'N'): 'peptide_bond', ('MET', 'CA', 'N'): 'peptide_bond', ('TRP', 'CH2', 'CZ2'): 'indole', ('ILE', 'CA', 'CB'): 'isoleucine', ('LYS', 'CD', 'CG'): 'lys_linker', ('PRO', 'CD', 'CG'): 'proline', ('GLU', 'C', 'O'): 'peptide_bond', ('TRP', 'CH2', 'CZ3'): 'indole', ('TPO', 'O1P', 'P'): 'phospho_thr', ('ILE', 'CA', 'N'): 'peptide_bond', ('CYS', 'C', 'CA'): 'peptide_bond', ('LYS', 'CD', 'CE'): 'lys_end', ('GLY', 'C', 'O'): 'peptide_bond', ('CYS', 'CA', 'CB'): 'cysteine', ('MET', 'C', 'CA'): 'peptide_bond', ('VAL', 'CA', 'N'): 'peptide_bond', ('THR', 'CB', 'CG2'): 'thr_linker', ('LEU', 'CA', 'CB'): 'leucine', ('ASP', 'C', 'O'): 'peptide_bond', ('LEU', 'CB', 'CG'): 'leucine', ('PRO', 'C', 'O'): 'peptide_bond', ('PHD', 'O1P', 'P'): 'phospho_asp', ('PHE', 'C', 'CA'): 'peptide_bond', ('ARG', 'CZ', 'NH2'): 'guanidinium', ('TYR', 'CB', 'CG'): 'tyr_linker', ('TYR', 'CE1', 'CZ'): 'phenyl', ('TYR', 'CE2', 'CZ'): 'phenyl', ('VAL', 'CB', 'CG2'): 'valine', ('LEU', 'C', 'O'): 'peptide_bond', ('LEU', 'CA', 'N'): 'peptide_bond', ('VAL', 'CA', 'CB'): 'valine', ('GLU', 'CA', 'N'): 'peptide_bond', ('LYS', 'C', 'O'): 'peptide_bond', ('SER', 'C', 'O'): 'peptide_bond', ('PTR', 'O2P', 'P'): 'phospho_tyr', ('GLU', 'CA', 'CB'): 'glu_linker', ('GLU', 'CB', 'CG'): 'glu_linker', ('GLN', 'CA', 'N'): 'peptide_bond', ('PHE', 'CD2', 'CG'): 'phenyl', ('TYR', 'C', 'CA'): 'peptide_bond', ('HIS', 'CA', 'CB'): 'his_linker', ('PTR', 'O1P', 'P'): 'phospho_tyr', ('TRP', 'CE2', 'CZ2'): 'indole', ('ASN', 'CG', 'OD1'): 'carboxamide', ('HIS', 'CA', 'N'): 'peptide_bond', ('ILE', 'CB', 'CG2'): 'isoleucine', ('HIS', 'CD2', 'NE2'): 'imidazole', ('TRP', 'CD2', 'CE3'): 'indole', ('ALA', 'C', 'CA'): 'peptide_bond', ('LEU', 'CD2', 'CG'): 'leucine', ('TRP', 'C', 'CA'): 'peptide_bond', ('GLN', 'CA', 'CB'): 'gln_linker', ('GLN', 'CB', 'CG'): 'gln_linker', ('HIP', 'OG', 'P'): 'phospho_his', ('THR', 'C', 'O'): 'peptide_bond', ('HIS', 'C', 'O'): 'peptide_bond', ('TRP', 'CD1', 'CG'): 'indole', ('ILE', 'C', 'CA'): 'peptide_bond', ('SER', 'CA', 'N'): 'peptide_bond', ('GLU', 'C', 'CA'): 'peptide_bond', ('GLN', 'C', 'CA'): 'peptide_bond', ('SEP', 'O1P', 'O3P'): 'phospho_ser', ('CYS', 'C', 'O'): 'peptide_bond', ('SER', 'CB', 'OG'): 'hydroxyl', ('HIS', 'CE1', 'NE2'): 'imidazole', ('PHE', 'CA', 'N'): 'peptide_bond', ('HIP', 'O1P', 'P'): 'phospho_his', ('ASP', 'CG', 'OD2'): 'carboxyl', ('ASN', 'CB', 'CG'): 'asn_linker', ('LYS', 'CE', 'NZ'): 'lys_end', ('CYS', 'CB', 'SG'): 'cysteine', ('VAL', 'CB', 'CG1'): 'valine', ('TYR', 'C', 'O'): 'peptide_bond', ('TPO', 'OG', 'P'): 'phospho_thr', ('VAL', 'C', 'O'): 'peptide_bond', ('HIP', 'O1P', 'O3P'): 'phospho_his', ('ARG', 'CD', 'NE'): 'guanidinium', ('ALA', 'C', 'O'): 'peptide_bond', ('LYS', 'CB', 'CG'): 'lys_linker', ('GLU', 'CD', 'OE2'): 'carboxyl', ('TPO', 'O3P', 'P'): 'phospho_thr', ('PHD', 'O2P', 'P'): 'phospho_asp', ('PHE', 'CA', 'CB'): 'phe_linker', ('ASN', 'CG', 'ND2'): 'carboxamide', ('HIS', 'CD2', 'CG'): 'imidazole', ('PHE', 'CD2', 'CE2'): 'phenyl', ('PRO', 'CA', 'N'): 'peptide_bond', ('ASN', 'C', 'O'): 'peptide_bond', ('ARG', 'CD', 'CG'): 'arg_linker', ('ARG', 'CA', 'N'): 'peptide_bond', ('TRP', 'CB', 'CG'): 'trp_linker', ('PHE', 'CB', 'CG'): 'phe_linker', ('HIS', 'C', 'CA'): 'peptide_bond', ('PHD', 'OG', 'P'): 'phospho_asp', ('MET', 'CE', 'SD'): 'met_end', ('GLN', 'C', 'O'): 'peptide_bond', ('HIS', 'CE1', 'ND1'): 'imidazole', ('CYS', 'CA', 'N'): 'peptide_bond', ('TRP', 'CE3', 'CZ3'): 'indole', ('TYR', 'CD2', 'CE2'): 'phenyl', ('ILE', 'CD1', 'CG1'): 'isoleucine', ('PTR', 'OG', 'P'): 'phospho_tyr', ('TRP', 'C', 'O'): 'peptide_bond', ('LYS', 'CA', 'N'): 'peptide_bond', ('TYR', 'CD1', 'CG'): 'phenyl', ('GLN', 'CD', 'CG'): 'gln_linker', ('ARG', 'CA', 'CB'): 'arg_linker', ('ASN', 'CA', 'N'): 'peptide_bond', ('SER', 'CA', 'CB'): 'ser_linker', ('GLN', 'CD', 'OE1'): 'carboxamide', ('PHD', 'O1P', 'O3P'): 'phospho_asp', ('THR', 'C', 'CA'): 'peptide_bond', ('ASN', 'CA', 'CB'): 'asn_linker', ('ASP', 'CB', 'CG'): 'asp_linker'}

SELECTION_NAME = 'query'  # The name for the query selection

# The following suffixes namespace the auto-generated objects.
# The Suns client reserves the right to freely delete anything unexpected within
# the following namespaces so that automated PyMOL clients can safely use these
# namespaces.
RESULT_SUFFIX  = 'result' # The suffix for each search result
FETCH_SUFFIX   = 'full'   # The suffix for each fetched context
SAVE_SUFFIX    = 'save'   # The suffix for saved results

SUNS_SERVER_ADDRESS = 'suns.degradolab.org' # The default message queue host

class SearchThread(threading.Thread):
    def __init__(
            self, rmsd, num_struct, random_seed, pdbstrs, server_address, cmd ):
        threading.Thread.__init__(self)
        self.request = json.dumps({
            'rmsd_cutoff'   : rmsd,
            'num_structures': num_struct,
            'random_seed'   : random_seed,
            'atoms'         : pdbstrs
        })
        # All PyMOL routines keep their own copy of the 'cmd' object.
        # I don't know why, but I copy that practice as a defensive precaution,
        # assuming they have a good reason for doing that.
        self.cmd = cmd
        
        self.channel = None # The AMQP Channel, unique per thread
        self.pdbs = {} # Tracks results to avoid duplicate generated names
        self.suns_server_address = server_address # Currently selected host
        
        self.begun = False # True once the thread begins consuming
        self.ended = False # True once the thread finishes consuming
         # This guards the 'begun' and 'end' variables
        self.lock = threading.Lock()
    
    def handle_delivery(
            self, channel, method_frame, header_frame, body, corr_id = None ):
        # Ignore messages with the wrong correlation ID so that we don't
        # contaminate the current results with results from an aborted previous
        # query
        if(corr_id == header_frame.correlation_id):
            # The inbound results don't use JSON, since there are a lot of them
            # and they tend to be small.
            
            # Message format is a tagged union and the first character is the
            # tag:
            # 0 = No more search results:      Payload = Empty
            # 1 = Search result:               Payload = PDB ID + PDB structure
            # 2 = Search time limit exceeded:  Payload = Empty
            # 3 = Error:                       Payload = message
            tag = body[0]
            if(tag == '0'):
                print '[*] Search done'
                self.stop()
            elif(tag == '1'):
                # The first 4 bytes of the payload are the PDB ID code
                pdbid = body[1:5]
                if(pdbid not in self.pdbs):
                    self.pdbs[pdbid] = 0
                sele_name = (
                    pdbid + '_%04d_%s' % (self.pdbs[pdbid], RESULT_SUFFIX))
                
                # If there is an existing result of the same name, delete it,
                # otherwise PyMOL will load the new result as an additional
                # model to the existing result.
                #
                # In principle there should never be an existing selection with
                # the same name since the client either migrates old results to
                # the <SAVE_SUFFIX> namespace or deletes unsaved ones when
                # launching a new search.  However, this is just a precaution in
                # case something goes wrong.
                self.cmd.delete(sele_name)
                
                # Load the structure into pymol.
                self.cmd.read_pdbstr(body[5:], sele_name)
                self.pdbs[pdbid] += 1
            elif(tag == '2'):
                print '[*] Time limit exceeded'
                self.stop()
            elif(tag == '3'):
                print '[*] Error: ' + body[1:]
            else:
                print '[*] Error: Invalid Server Response'
                print '[*] Response = \'' + body + '\''
    
    # Send a search request and await results
    def run(self):
        # This correlation ID ensures that each thread does not receive results
        # from previous search requests that prematurely aborted
        corr_id = str(uuid.uuid4())
        
        credentials = PlainCredentials('suns-client', 'suns-client')
        try:
            # Connection initialization
            connection = BlockingConnection(
                ConnectionParameters(
                    host = self.suns_server_address,
                    credentials = credentials,
                    virtual_host = 'suns-vhost'))
            try:
                # Channel Initialization
                self.channel = connection.channel()
                self.channel.exchange_declare(
                    exchange = 'suns-exchange-responses',
                    passive = True,
                    durable = True)
                self.channel.exchange_declare(
                    exchange = 'suns-exchange-requests',
                    passive = True,
                    durable = True)
                
                # Queue Initialization
                result = self.channel.queue_declare(exclusive = True)
                try:
                    # Auto-generates an anonymous queue name
                    callback_queue = result.method.queue
                    self.channel.queue_bind(
                        exchange = 'suns-exchange-responses',
                        queue = callback_queue,
                        routing_key = callback_queue)
                    self.channel.basic_consume(
                        lambda c, m, h, b :
                            self.handle_delivery(c, m, h, b, corr_id),
                        no_ack = True,
                        queue = callback_queue)
                    
                    # Send Request
                    print '[*] Waiting for Server...'
                    self.channel.basic_publish(
                        exchange = 'suns-exchange-requests',
                        routing_key = '1.0.0',
                        properties = BasicProperties(
                            reply_to = callback_queue,
                            correlation_id = corr_id),
                        body = self.request)
                    
                    try:
                        self.lock.acquire()
                        try:
                            self.begun = True
                            ended = self.ended
                        finally:
                            self.lock.release()
                        # Race condition:
                        #     This thread could receive a stop() signal before
                        #     start_consuming or could receive a second
                        #     stop() signal after done consuming.  Pika does not
                        #     expose a good way to release the lock immediately
                        #     after beginning consuming or to acquire the lock
                        #     immediately before done consuming.
                        #
                        # The chance of triggering this race condition is very
                        # low and the user can fix it by restarting the wizard
                        if (not ended):
                            # This call is a blocking call.
                            # It blocks until channel.stop_consuming is called.
                            self.channel.start_consuming()
                    finally:
                        self.lock.acquire()
                        try:
                            self.ended = True
                        finally:
                            self.lock.release()
                        self.cmd.orient(SELECTION_NAME)
                finally:
                    self.channel.queue_delete(queue = callback_queue)
            finally:
                connection.close()
        except socket.error:
            print ("[*] Error: Unable to connect to '"
                 + self.suns_server_address
                 + "'")
    
    def cancel_search(self):
        print '[*] Search cancelled'
        self.stop()
    
    def stop(self):
        self.lock.acquire()
        try:
            if(not self.ended):
                if(self.begun): # The channel should not be None if begun = True
                    self.channel.stop_consuming()
                self.ended = True
        finally:
            self.lock.release()


# This wizard requires the ability to select bonds instead of atoms, since all
# motifs are defined by bonds.  The only way I know how to do that is to
# override Python's "PkTB" input.  I document various vagaries of that in the
# appropriate sections below.
class Suns_search(Wizard):
    def __init__(self, _self = cmd):
        Wizard.__init__(self, _self)
        
        # Clean the slate
        self.cmd.unpick()        
        self.do_select(SELECTION_NAME, 'none')
        self.prompt = ['Select motifs using left click']
        
        # type Key = (Object, Model, Segi, Chain, Resn, Resi, Name)
        # type Val = Selection corresponding to the key
        # word_list :: Dict Key Val
        self.word_list = {}
        
        # Default values
        self.rmsd_cutoff = 1.0
        self.random_seed = 0
        self.number_of_structures = 100
        self.suns_server_address = SUNS_SERVER_ADDRESS
        
        self.searchThread = None
        
        # Rebind left click to select bonds.  I choose left click since it is
        # the most portable input method across devices, especially laptops, and
        # also because the user definitely does not need the old left click
        # function for the duration of the wizard.
        self.cmd.button('single_left', 'None', 'PkTB')
        
        # PyMOL has this really obscure and undocumented 'auto_hide_selections'
        # setting that auto hides selections when using editing commands.  This
        # wizard overrides "PkTB", which is one of those editing commands, but
        # we still want to display selections to the user, so we must disable
        # this setting.
        self.prev_auto_hide_setting = self.cmd.get('auto_hide_selections')
        self.cmd.set('auto_hide_selections',0)
    
    def cleanup(self):
        self.stop_search()
        
        # This is a hack.  I don't yet know how to retrieve the previous value
        # of this button setting so that I can restore the original value.  For
        # now I assume that the user was in a viewing mode, where the default
        # behavior for left click is '+/-'.
        #
        # Very few users use anything other than viewing mode, especially in
        # conjunction with searching.  Also, users can fix the mistaken
        # reversion by cycling through mouse button modes, which will properly
        # reset the action for single left click.
        self.cmd.button('single_left', 'None', '+/-')
        
        self.cmd.set('auto_hide_selections', self.prev_auto_hide_setting)
        self.cmd.delete(SELECTION_NAME)
    
    def set_rmsd(self, rmsd):
        self.rmsd_cutoff = rmsd
        self.cmd.refresh_wizard()
    
    def set_num_structures(self, num_structures):
        self.number_of_structures = num_structures
        self.cmd.refresh_wizard()
    
    def ask_random_seed(self, ask_for_random_seed):
        if(ask_for_random_seed):
            try:
                self.random_seed = int(tkSimpleDialog.askstring(
                    'Randomize Order','Random Seed (Integer):'))
            except:
                self.random_seed = 0
        else:
            self.random_seed = 0
        self.cmd.refresh_wizard()
    
    def create_random_seed_menu(self):
        random_seed_menu = [[2, 'Randomize Order', '']]
        random_seed_menu.append(
            [1, 'No', 'cmd.get_wizard().ask_random_seed(False)'])
        random_seed_menu.append(
            [1, 'Yes', 'cmd.get_wizard().ask_random_seed(True)'])
        return random_seed_menu
        
    def create_rmsd_menu(self):
        rmsd_menu = [[2, 'RMSD Cutoff', '']]
        for rmsd_choice in range(1,21):
            rmsd = float(rmsd_choice) / 10.0
            rmsd_menu.append(
                [1, str(rmsd) , 'cmd.get_wizard().set_rmsd(' + str(rmsd) + ')'])
        return rmsd_menu
    
    def create_num_structures_menu(self):
        num_structures_menu = [[2, 'Number of Results', '']]
        for n in [10, 20, 50, 100, 200, 500, 1000]:
            num_structures_menu.append(
                [1, str(n), 'cmd.get_wizard().set_num_structures(' + str(n) + ')'])
        return num_structures_menu
    
    def create_server_address_menu(self):
        server_address_menu = [[2, 'Server Address', '']]
        server_address_menu.append(
            [1, 'Default: ' + SUNS_SERVER_ADDRESS, 'cmd.get_wizard().set_server_address("' + SUNS_SERVER_ADDRESS + '")'])
        server_address_menu.append(
            [1, 'User defined', 'cmd.get_wizard().ask_server_address()'])
        return server_address_menu
        
    def set_server_address(self, server_address):
        self.suns_server_address = server_address
        self.cmd.refresh_wizard()
        
    def ask_server_address(self):
        try:
            server_address = tkSimpleDialog.askstring(
                'Suns server address','Suns server address:')
        except:
            server_address = SUNS_SERVER_ADDRESS
        self.set_server_address(server_address)
    
    def get_panel(self):
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
            [ 2, 'Clear Results', 'cmd.get_wizard().delete_results()'],
            [ 2, 'Clear Selection','cmd.get_wizard().clear_selection()'],
            [ 2, 'Clear Saved', 'cmd.get_wizard().delete_saved()'],
            [ 2, 'Fetch Full Contexts','cmd.get_wizard().fetch_full_context()'],
            [ 2, 'Clear Contexts','cmd.get_wizard().delete_full()'],
            [ 2, 'Done','cmd.set_wizard()'] ]
    
    def fetch_full_context(self):
        '''
        This method will look at what search results are currently enabled,
        and then fetch and align the full structure so that the user can 
        see the full context of the match.
        '''
        objs = self.cmd.get_names('objects', enabled_only = 1)
        
        def match(objName):
            words = objName.split('_')
            return (len(words) == 3) and (words[2] == RESULT_SUFFIX)
        
        matches = filter(match, objs)
        if(len(matches) > 1):
            if(not tkMessageBox.askyesno(
                "Multiple Results Selected",
                "Fetching multiple contexts takes time and will block PyMOL until finished.  Proceed?")):
                return;
        
        for obj in matches:
            w = obj.split('_')
            del w[-1]
            w.append(FETCH_SUFFIX)
            new_object_name = '_'.join(w)
            self.cmd.fetch(w[0], new_object_name)
            attr_dict = {'x' : []}
            # Get the atom info for the entire result, not just the part that
            # matched the query, in order to improve the alignment
            self.cmd.iterate(
                obj,
                'x.append( [model,segi,chain,resn,resi,name,alt] )',
                space=attr_dict)
            selection = []
            for item in attr_dict['x']:
                altLoc = item[6]
                # PyMOL chokes on empty arguments to 'alt'
                if(altLoc.strip() == ''):
                    selection += [
                        '(chain %s and resn %s and resi %s and name %s)' % tuple(item[2:6])]
                else:
                    selection += [
                        '(chain %s and resn %s and resi %s and name %s and alt %s)' % tuple(item[2:7])]
            selection = '(' + ' or '.join(selection) + ')'
            self.cmd.pair_fit(
                new_object_name + ' and (' + selection + ')',
                obj             + ' and (' + selection + ')')
        self.cmd.orient(SELECTION_NAME)
    
    def delete_results(self):
        self.cmd.delete('*_' + RESULT_SUFFIX)
    
    def delete_saved(self):
        self.cmd.delete('*_' + SAVE_SUFFIX)
    
    def delete_full(self):
        self.cmd.delete('*_' + FETCH_SUFFIX)
    
    def launch_search(self):
        pdbstr = self.cmd.get_pdbstr(SELECTION_NAME)
        self.stop_search()
        self.delete_results()
        self.searchThread = SearchThread(self.rmsd_cutoff, self.number_of_structures, self.random_seed, pdbstr, self.suns_server_address, self.cmd)
        self.searchThread.start()
    
    def clear_selection(self):
        # Do not delete the selection, otherwise PyMOL will throw an exception
        # when you try to search.  Instead, set it to the empty selection.
        self.word_list = {}
        self.do_select(SELECTION_NAME, 'none')
    
    def cancel_search(self):
        if(self.searchThread is not None):
            self.searchThread.cancel_search()
    
    def stop_search(self):
        if(self.searchThread is not None):
            self.searchThread.stop()
    
    def do_select(self, name, selection_logic):
        if(selection_logic != ''):
            self.cmd.select(name, selection_logic)
            self.cmd.enable(name)
        else:
            self.cmd.select(name, 'none')
            self.cmd.disable(name)
            
    def do_pick(self, bondFlag):
        '''
        This is the method that is called each time the user
        uses the mouse to select a bond or atom.
        '''
        # I think bondFlag only = 1 if we are selecting bonds.
        # We are only accepting bond selections.
        if(bondFlag == 1):
            attr_dict = {'x' : []}
            # Get the atom info for the two atoms of the bond.
            obj, unused = self.cmd.index("pk1")[0]
            words = obj.split('_')
            if(len(words) == 3 and words[2] == RESULT_SUFFIX):
                del words[2]
                words.append(SAVE_SUFFIX)
                newObj = '_'.join(words)
                self.cmd.set_name(obj, newObj)
                obj = newObj
            self.cmd.iterate("pk1", 'x.append( (model,segi,chain,resn,resi,name,alt) )', space = attr_dict)
            self.cmd.iterate("pk2", 'x.append( (model,segi,chain,resn,resi,name,alt) )', space = attr_dict)
            bond = [attr_dict['x'][0][5], attr_dict['x'][1][5]]
            bond.sort()
            
            # We only handle bonds per residue, no residue spanning bonds.
            if( (attr_dict['x'][0][2] != attr_dict['x'][1][2]) or (attr_dict['x'][0][4] != attr_dict['x'][1][4]) ):
                return
            
            # The key is of the form (resn, atom0, atom1) where the atom names are in alphabetical order.
            key = (attr_dict['x'][0][3], bond[0], bond[1])
            # Look up what word this bond is part of.
            if(key in BOND_WORD_DICT):
                word = BOND_WORD_DICT[key]
                # Now form a new key that has the current residue and the word.
                key = tuple([obj] + list(attr_dict['x'][0][0:5]) + [word])
                if(key in self.word_list):
                    # Disable the selection
                    del self.word_list[key]
                else:
                    # Enable the selection
                    if(key[2].strip() == ''):
                        self.word_list[key] = '(object %s and model %s and chain %s and resn %s and resi %s and name %s )' % tuple(list(key[0:2]) + list(key[3:6]) + [WORDS_DICT[word][key[4]]])
                    else:
                        self.word_list[key] = '(object %s and model %s and segi %s and chain %s and resn %s and resi %s and name %s )' % tuple(list(key[0:6]) + [WORDS_DICT[word][key[4]]])
            #else: # If this isn't part of a word, we'll still let them select the bond.
            #    # If this bond isn't part of a word, then the key
            #    # will include the full identifier for the two atoms.
            #    # We should sort the atoms so that if the bond comes in as a,b or b,a
            #    # we'll know it's the same thing.  For now, I'll just sort on the atom name.
            #    # This could cause a problem if the two atoms have the same name, but really,
            #    # should a bond have two atoms of the same name?
            #    if(attr_dict['x'][0][5] < attr_dict['x'][1][5]):
            #        key = tuple(list(attr_dict['x'][0][0:6]) + list(attr_dict['x'][1][0:6]))
            #    else:
            #        key = tuple(list(attr_dict['x'][1][0:6]) + list(attr_dict['x'][0][0:6]))
            #    if(key in self.word_list):
            #        del self.word_list[key]
            #    else:
            #        self.word_list[key] = '((model %s and segi %s and chain %s and resn %s and resi %s and name %s) or (model %s and segi %s and chain %s and resn %s and resi %s and name %s) )' % key
        
        self.cmd.unpick()
        self.do_select(SELECTION_NAME, ' or '.join(self.word_list.values()).strip())
