import re
from Bio import SeqIO
import gzip
from Bio.Seq import Seq
from inimotif_core import KmerCounter
from yattag import Doc,indent
from windows import gen_full_win_list
import warnings


# a pattern for matching DNA repeats
class RepeatPattern:
    def __init__(self, seq, n_min_rep, revcom_flag=True, n_max_rep=None):
        self.revcom_flag = revcom_flag
        
        if n_max_rep is None:
            n_max_rep = ''
        self.forward_pattern = re.compile(f'({seq}){{{n_min_rep},{n_max_rep}}}')
        
        if self.revcom_flag:
            kmer_counter = KmerCounter(len(seq))
            revcom_seq = kmer_counter.revcom(seq)
            self.revcom_pattern = re.compile(f'({revcom_seq}){{{n_min_rep},{n_max_rep}}}')
    
    def mask(self, in_str):
        def myrepl(m):
            return ('N'*len(m.group(0)))
        in_str = re.sub(self.forward_pattern, myrepl, in_str)
        if self.revcom_flag:
            in_str = re.sub(self.revcom_pattern, myrepl, in_str)
        return in_str

class Motif:
    def __init__(self, seq, n_max_mutation=0, revcom_flag=True):
        assert n_max_mutation<len(seq), f'n_max_mutation={n_max_mutation} is smaller than seq length {len(seq)}!'
        kc = KmerCounter(len(seq))
        self.kc = kc
        self.revcom_flag = revcom_flag
        
        self.seq = seq
        self.n_max_mutation = n_max_mutation
        
        forward_seq_hash = kc.kmer2hash(seq)
        self.forward_hamball = kc.get_hamming_ball(forward_seq_hash,n_max_mutation)
        
        self.revcom_flag = revcom_flag
        if revcom_flag:
            revcom_seq_hash = kc.revcom_hash(forward_seq_hash)
            self.is_palindrome = forward_seq_hash==revcom_seq_hash
            self.revcom_seq = kc.hash2kmer(revcom_seq_hash)
            self.revcom_hamball = kc.get_hamming_ball(revcom_seq_hash,n_max_mutation)
    
    def __str__(self):
        if not self.revcom_flag:
            return f'{self.seq} (Forward) n_max_mutation={self.n_max_mutation}'
        elif self.is_palindrome:
            return f'{self.seq} (Palindrome) n_max_mutation={self.n_max_mutation}'
        else:
            return f'{self.seq} (Forward) {self.revcom_seq} (Revcom) n_max_mutation={self.n_max_mutation}'
    
    def mask(self, in_str):
        in_str_list = [c for c in in_str]
        len_str = len(in_str)
        k = self.kc.k

        i=0
        prev_hash = self.kc.dtype(-1)
        while(i<=len_str-k):
            tmpstr = in_str[i:i+k]
            # omit a kmer if it contains "N"
            if "N" in tmpstr:
                i += 1 + max([pos for pos,char in enumerate(tmpstr) if char == 'N'])
                prev_hash = self.kc.dtype(-1)
                continue
            if prev_hash==self.kc.dtype(-1):
                tmphash = self.kc.kmer2hash(tmpstr)
            # reuse hash in previous position
            else:
                tmphash = (prev_hash<< self.kc.dtype(2) ) & self.kc.mask
                tmphash += self.kc.base[ in_str[i+k-1] ]
            prev_hash = tmphash
            forward_flag = tmphash in self.forward_hamball
            revcom_flag = self.revcom_flag and (tmphash in self.revcom_hamball)
            if forward_flag or revcom_flag:
                for ii in range(i,i+k):
                    in_str_list[ii] = 'N' 
                i += k
                continue
            i += 1
        return "".join(in_str_list)
    
    # scan motif in input string and report its locations
    def scan(self, in_str):
        pos_list = []
        len_str = len(in_str)
        k = self.kc.k

        i=0
        prev_hash = self.kc.dtype(-1)
        while(i<=len_str-k):
            tmpstr = in_str[i:i+k]
            # omit a kmer if it contains "N"
            if "N" in tmpstr:
                i += 1 + max([pos for pos,char in enumerate(tmpstr) if char == 'N'])
                prev_hash = self.kc.dtype(-1)
                continue
            if prev_hash==self.kc.dtype(-1):
                tmphash = self.kc.kmer2hash(tmpstr)
            # reuse hash in previous position
            else:
                tmphash = (prev_hash<< self.kc.dtype(2) ) & self.kc.mask
                tmphash += self.kc.base[ in_str[i+k-1] ]
            prev_hash = tmphash
            forward_flag = tmphash in self.forward_hamball
            revcom_flag = self.revcom_flag and (tmphash in self.revcom_hamball)
            if forward_flag or revcom_flag:
                pos_list.append(i)
                i += k
                continue
            i += 1
        return pos_list

class Masker:
    def __init__(self):
        self.pattern_list = []
    
    def clear(self):
        self.pattern_list = []
    
    # add repetitive pattern    
    def add_reppat(self, seq, n_min_rep, revcom_flag):
        self.pattern_list.append( RepeatPattern(seq, n_min_rep, revcom_flag) )
    
    # add motif
    def add_motif(self, seq, n_max_mutation, revcom_flag):
        self.pattern_list.append( Motif(seq, n_max_mutation, revcom_flag) )
        
    def mask(self,in_str):
        in_str = in_str.upper()
        for pat in self.pattern_list:
            in_str = pat.mask(in_str)
        return in_str
    
    def mask_file(self,input_fasta_file_name, out_file=None):
        if out_file is None:
            out_file = f'{input_fasta_file_name}.mask.fasta'
            
        foh = open(out_file,'w')
        if input_fasta_file_name.endswith(".gz"):
            fh = gzip.open(input_fasta_file_name,"rt")
        else:
            fh = open(input_fasta_file_name, "r")
            
        for rec in SeqIO.parse(fh,"fasta"):
            tmpstr = self.mask(str(rec.seq))
            rec.seq = Seq(tmpstr)
            SeqIO.write(rec,foh,'fasta')
        
        foh.close()
        fh.close()
        
        

class MotifScanner:
    def __init__(self):
        self.motif_list = []
    
    def clear(self):
        self.motif_list = []
    
    def add_motif(self, seq, n_max_mutation, revcom_flag):
        self.motif_list.append( Motif(seq, n_max_mutation, revcom_flag) )
        
    def scan(self, in_str):
        res_list = [m.scan(in_str) for m in self.motif_list]  # nested list
        return res_list
    
    def _get_style_str(self):
        style_str = """
        body {
            font-family: sans-serif, Arial, Helvetica;
        }
        p {
            word-break: break-all;
            white-space: normal;
        }
        motif1{
            color: red;
        }
        motif2{
            color: blue;
        }
        overlap{
            color: green;
        }
        """
        return style_str
    
    def scan_file(self, input_fasta_file_name, out_file="motif_scan.html"):
        if len(self.motif_list)>2:
            self.motif_list = self.motif_list[0:2]
            warnings.warn('Currently only support scanning first two motifs')
        if not out_file.endswith(".html"):
            out_file += ".html"
            warnings.warn(f'add html suffix to out_file="{out_file}"')
        
        if input_fasta_file_name.endswith(".gz"):
            fh = gzip.open(input_fasta_file_name,"rt")
        else:
            fh = open(input_fasta_file_name, "r")
        style_str = self._get_style_str()
        doc, tag, text = Doc().tagtext()
        doc.asis('<!DOCTYPE html>')
        with tag('html',lang="en"):
            with tag('head'):
                with tag('title'):
                    text('Motif Sequences')
                doc.stag('meta', charset="utf-8")
                doc.stag('meta', name="viewport", content="width=device-width, initial-scale=1")
                with tag('style'):
                    text(style_str) # specify style
            with tag('body'):
                with tag('h1'):
                    text('IniMotif')
                with tag('h2'):
                    text('Sequences contain motifs')
                    doc.stag('br')
                    doc.stag('br')
                    for i,motif in enumerate(self.motif_list):
                        with tag(f'motif{i+1}'):
                            text(f'motif {i+1}: {str(motif)}')
                            doc.stag('br')
                for rec in SeqIO.parse(fh,"fasta"):
                    tmpseq = str(rec.seq).upper()
                    m1_pos_list,m2_pos_list = self.scan(tmpseq)
                    if not m1_pos_list and not m2_pos_list:
                        continue
                    
                    win_list,win_type_list = gen_full_win_list(m1_pos_list, m2_pos_list,
                                                                len(self.motif_list[0].seq),len(self.motif_list[1].seq), len(tmpseq))

                    with tag('p'):
                        # output header
                        text(">"+rec.id)
                        # output line break
                        doc.stag('br')
                        # output sequence
                        for win,win_type in zip(win_list, win_type_list):
                            if win_type==0:
                                text(tmpseq[win[0]:win[1]])
                            elif win_type==1:
                                with tag('motif1'):
                                    text(tmpseq[win[0]:win[1]])
                            elif win_type==2:
                                with tag('motif2'):
                                    text(tmpseq[win[0]:win[1]])
                            elif win_type==3:
                                with tag('overlap'):
                                    text(tmpseq[win[0]:win[1]])
                            else:
                                warnings.warn(f'Unkown win_type={win_type}')
        fh.close()

        html_str = indent(doc.getvalue(), indent_text = True) # will also indent the text directly contained between <tag> and </tag>
        with open(out_file,'w') as out_fh:
            out_fh.write(html_str)

    
if __name__=="__main__":
    in_str = 'AAAAAAAAAAACGTGCCCCCGTGGGGGCGTGAAAACACGCCCCCACGTTTTCACGTTTTCACG'
    mask = Masker()
    mask.add_reppat('A',4,True)
    mask.add_reppat('C',4,True)
    mask.add_motif('CGTG',1,False)
    # mask.add_motif_pat('CGTG',1,True)
    
    # out_str = mask.mask(in_str)
    # print(in_str)
    # print(out_str)
    
    infile = './exampledata/NF1-1'
    outfile = './test_mask.fasta'
    mask.mask_file(infile, outfile)
    
    ms = MotifScanner()
    ms.add_motif('CGTG',1,True)
    ms.add_motif('AAAA',1,False)
    # print([(i,c) for i,c in enumerate(in_str)])
    # print(ms.scan(in_str))
    
    nfile = './exampledata/NF1-1'
    out_file = './test_motif_scan.html'
    ms.scan_file(infile, out_file)
    
    