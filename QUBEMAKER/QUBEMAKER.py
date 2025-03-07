# class to read the protein pdb files and store the infomation?

from collections import OrderedDict
from re import sub



class QUBEMAKER:
    """"A class to make the required QUBE-OPEMMM system for a protein ligand complex.
    The class contains templates for each residue that are combined into a specific protein xml file.
    This is due to repeted residues having different non-bonded sections which must be assigned.
    Write the OPLS combination rule to the xml file?"""

    def __init__(self, PDB_file, PSF_file=None, PAR_file=None):
        self.PDB = PDB_file
        self.PSF = PSF_file
        self.PAR = PAR_file
        self.element_dict = {'H': 1.008000,  # Group 1
                        'C': 12.011000,  # Group 4
                        'N': 14.007000, 'P': 30.973762,  # Group 5
                        'O': 15.999000, 'S': 32.060000,  # Group 6
                        'F': 18.998403, 'Cl': 35.450000, 'Br': 79.904000, 'I': 126.904470  # Group 7
                        }
        self.natom = 0
        self.molecule = []
        self.charges = []
        self.esp = []
        self.sigma = []
        self.atom_names = []
        self.read_pdb()
        self.write_parameters()


    def write_parameters(self):
        """Take the molecule's parameter set and write an xml file for the molecule."""

        # write the standard parameter file
        self.write_QUBE_xml()

        # now get the esp and sigma values in the order of the pdb
        self.read_par()

        # now get the charges in the order of the pdb
        # self.read_psf()

        # now write to file in the script

    def read_par(self):
        """Read the par file to get the L-J parameters."""

        with open(self.PAR, 'r') as par:
            lines = par.readlines()

        for i, line in enumerate(lines):
            if 'NONBONDED' in line:
                for line in lines[i+9:]:
                    print(line)
                break


    def read_pdb(self):

        with open(self.PDB, 'r') as pdb:
            lines = pdb.readlines()

        molecule = []

        # atom counter used for graph node generation
        for line in lines:
            if 'ATOM' in line or 'HETATM' in line:
                element = str(line[76:78])
                element = sub('[0-9]+', '', element)
                element = element.replace(" ", "")
                self.atom_names.append(str(line.split()[2]))

                # If the element column is missing from the pdb, extract the element from the name.
                if not element:
                    element = str(line.split()[2])[:-1]
                    element = sub('[0-9]+', '', element)

                # Also add the atom number as the node in the graph
                self.natom += 1
                molecule.append([element, float(line[30:38]), float(line[38:46]), float(line[46:54])])
        self.molecule = molecule
        print(self.molecule)

    # def read_pdb(self):
    #     """Read the pdb file and collect the residue infomation."""
    #
    #
    #     with open(self.PDB, 'r') as pdb:
    #         lines = pdb.readlines()
    #
    #     # start parsing through the lines of the pdb and extract the infomation
    #     for line in lines:
    #         if 'ATOM' in line or 'HETATM' in line:
    #             # found out the element of the atom
    #             element = str(line[76:78])
    #             element = sub('[0-9]+', '', element)
    #             element = element.replace(" ", "")
    #             print(element)
    #
    #             # get the residue name
    #             residue = str(line[17:20])
    #             print(residue)
    #
    #             # get the residue number
    #             residue_no = str(line[25:27])
    #             print(residue_no)
    #
    #             # get the atom name
    #             atom_name = str(line.split()[2])
    #
    #             # store the name as index, element and the mass of the atom other info added later
    #             self.atoms[self.natom] = [element, self.element_dict[element], atom_name]
    #
    #             # store the residue infomation
    #             if residue_no in self.residue.keys():
    #                 self.residue[residue_no][1].append(atom_name)
    #             else:
    #                  self.residue[residue_no] = [residue, [atom_name]]
    #
    #
    #             # add one for the next atom type
    #             self.natom += 1
    #
    #     print(f'PDB read and {len(self.residue)} residues found with sequence:')
    #     for key in self.residue.keys():
    #         print(f'{self.residue[key][0]}')
    #     print(self.residue)
    #     print(self.atoms)

    def write_QUBE_xml(self):
        """This function builds the standard QUBE.xml which is then edited to add the custom nonbonded params."""
        with open('QUBEMAKER.xml', 'w+')as xml:
            xml.write('''<ForceField>
 <AtomTypes>
  <Type name="0" class="N000" element="N" mass="14.00672"/>
  <Type name="1" class="H001" element="H" mass="1.007947"/>
  <Type name="2" class="C002" element="C" mass="12.01078"/>
  <Type name="3" class="H002" element="H" mass="1.007947"/>
  <Type name="4" class="C004" element="C" mass="12.01078"/>
  <Type name="5" class="H005" element="H" mass="1.007947"/>
  <Type name="6" class="C006" element="C" mass="12.01078"/>
  <Type name="7" class="O007" element="O" mass="15.99943"/>
  <Type name="8" class="N008" element="N" mass="14.00672"/>
  <Type name="9" class="H009" element="H" mass="1.007947"/>
  <Type name="10" class="C010" element="C" mass="12.01078"/>
  <Type name="11" class="H011" element="H" mass="1.007947"/>
  <Type name="12" class="C012" element="C" mass="12.01078"/>
  <Type name="13" class="H013" element="H" mass="1.007947"/>
  <Type name="14" class="C014" element="C" mass="12.01078"/>
  <Type name="15" class="H015" element="H" mass="1.007947"/>
  <Type name="16" class="C016" element="C" mass="12.01078"/>
  <Type name="17" class="H017" element="H" mass="1.007947"/>
  <Type name="18" class="N018" element="N" mass="14.00672"/>
  <Type name="19" class="H019" element="H" mass="1.007947"/>
  <Type name="20" class="C020" element="C" mass="12.01078"/>
  <Type name="21" class="N021" element="N" mass="14.00672"/>
  <Type name="22" class="H022" element="H" mass="1.007947"/>
  <Type name="23" class="N023" element="N" mass="14.00672"/>
  <Type name="24" class="H024" element="H" mass="1.007947"/>
  <Type name="25" class="C025" element="C" mass="12.01078"/>
  <Type name="26" class="O026" element="O" mass="15.99943"/>
  <Type name="27" class="N027" element="N" mass="14.00672"/>
  <Type name="28" class="H028" element="H" mass="1.007947"/>
  <Type name="29" class="C029" element="C" mass="12.01078"/>
  <Type name="30" class="H030" element="H" mass="1.007947"/>
  <Type name="31" class="C031" element="C" mass="12.01078"/>
  <Type name="32" class="H032" element="H" mass="1.007947"/>
  <Type name="33" class="C033" element="C" mass="12.01078"/>
  <Type name="34" class="O034" element="O" mass="15.99943"/>
  <Type name="35" class="O035" element="O" mass="15.99943"/>
  <Type name="36" class="H036" element="H" mass="1.007947"/>
  <Type name="37" class="C037" element="C" mass="12.01078"/>
  <Type name="38" class="O038" element="O" mass="15.99943"/>
  <Type name="39" class="N039" element="N" mass="14.00672"/>
  <Type name="40" class="H040" element="H" mass="1.007947"/>
  <Type name="41" class="C041" element="C" mass="12.01078"/>
  <Type name="42" class="H042" element="H" mass="1.007947"/>
  <Type name="43" class="C043" element="C" mass="12.01078"/>
  <Type name="44" class="H044" element="H" mass="1.007947"/>
  <Type name="45" class="C045" element="C" mass="12.01078"/>
  <Type name="46" class="O046" element="O" mass="15.99943"/>
  <Type name="47" class="N047" element="N" mass="14.00672"/>
  <Type name="48" class="H048" element="H" mass="1.007947"/>
  <Type name="49" class="C049" element="C" mass="12.01078"/>
  <Type name="50" class="O050" element="O" mass="15.99943"/>
  <Type name="51" class="N051" element="N" mass="14.00672"/>
  <Type name="52" class="H052" element="H" mass="1.007947"/>
  <Type name="53" class="C053" element="C" mass="12.01078"/>
  <Type name="54" class="H054" element="H" mass="1.007947"/>
  <Type name="55" class="C055" element="C" mass="12.01078"/>
  <Type name="56" class="H056" element="H" mass="1.007947"/>
  <Type name="57" class="C057" element="C" mass="12.01078"/>
  <Type name="58" class="O058" element="O" mass="15.99943"/>
  <Type name="59" class="O059" element="O" mass="15.99943"/>
  <Type name="60" class="C060" element="C" mass="12.01078"/>
  <Type name="61" class="O061" element="O" mass="15.99943"/>
  <Type name="62" class="N062" element="N" mass="14.00672"/>
  <Type name="63" class="H063" element="H" mass="1.007947"/>
  <Type name="64" class="C064" element="C" mass="12.01078"/>
  <Type name="65" class="H065" element="H" mass="1.007947"/>
  <Type name="66" class="C066" element="C" mass="12.01078"/>
  <Type name="67" class="H067" element="H" mass="1.007947"/>
  <Type name="68" class="S068" element="S" mass="32.0655"/>
  <Type name="69" class="C069" element="C" mass="12.01078"/>
  <Type name="70" class="O070" element="O" mass="15.99943"/>
  <Type name="71" class="N071" element="N" mass="14.00672"/>
  <Type name="72" class="H072" element="H" mass="1.007947"/>
  <Type name="73" class="C073" element="C" mass="12.01078"/>
  <Type name="74" class="H074" element="H" mass="1.007947"/>
  <Type name="75" class="C075" element="C" mass="12.01078"/>
  <Type name="76" class="H076" element="H" mass="1.007947"/>
  <Type name="77" class="S077" element="S" mass="32.0655"/>
  <Type name="78" class="H078" element="H" mass="1.007947"/>
  <Type name="79" class="C079" element="C" mass="12.01078"/>
  <Type name="80" class="O080" element="O" mass="15.99943"/>
  <Type name="81" class="N081" element="N" mass="14.00672"/>
  <Type name="82" class="H082" element="H" mass="1.007947"/>
  <Type name="83" class="C083" element="C" mass="12.01078"/>
  <Type name="84" class="H084" element="H" mass="1.007947"/>
  <Type name="85" class="C085" element="C" mass="12.01078"/>
  <Type name="86" class="H086" element="H" mass="1.007947"/>
  <Type name="87" class="S087" element="S" mass="32.0655"/>
  <Type name="88" class="C088" element="C" mass="12.01078"/>
  <Type name="89" class="O089" element="O" mass="15.99943"/>
  <Type name="90" class="N090" element="N" mass="14.00672"/>
  <Type name="91" class="H091" element="H" mass="1.007947"/>
  <Type name="92" class="C092" element="C" mass="12.01078"/>
  <Type name="93" class="H093" element="H" mass="1.007947"/>
  <Type name="94" class="C094" element="C" mass="12.01078"/>
  <Type name="95" class="H095" element="H" mass="1.007947"/>
  <Type name="96" class="C096" element="C" mass="12.01078"/>
  <Type name="97" class="H097" element="H" mass="1.007947"/>
  <Type name="98" class="C098" element="C" mass="12.01078"/>
  <Type name="99" class="O099" element="O" mass="15.99943"/>
  <Type name="100" class="O100" element="O" mass="15.99943"/>
  <Type name="101" class="H101" element="H" mass="1.007947"/>
  <Type name="102" class="C102" element="C" mass="12.01078"/>
  <Type name="103" class="O103" element="O" mass="15.99943"/>
  <Type name="104" class="N104" element="N" mass="14.00672"/>
  <Type name="105" class="H105" element="H" mass="1.007947"/>
  <Type name="106" class="C106" element="C" mass="12.01078"/>
  <Type name="107" class="H107" element="H" mass="1.007947"/>
  <Type name="108" class="C108" element="C" mass="12.01078"/>
  <Type name="109" class="H109" element="H" mass="1.007947"/>
  <Type name="110" class="C110" element="C" mass="12.01078"/>
  <Type name="111" class="H111" element="H" mass="1.007947"/>
  <Type name="112" class="C112" element="C" mass="12.01078"/>
  <Type name="113" class="O113" element="O" mass="15.99943"/>
  <Type name="114" class="N114" element="N" mass="14.00672"/>
  <Type name="115" class="H115" element="H" mass="1.007947"/>
  <Type name="116" class="C116" element="C" mass="12.01078"/>
  <Type name="117" class="O117" element="O" mass="15.99943"/>
  <Type name="118" class="N118" element="N" mass="14.00672"/>
  <Type name="119" class="H119" element="H" mass="1.007947"/>
  <Type name="120" class="C120" element="C" mass="12.01078"/>
  <Type name="121" class="H121" element="H" mass="1.007947"/>
  <Type name="122" class="C122" element="C" mass="12.01078"/>
  <Type name="123" class="H123" element="H" mass="1.007947"/>
  <Type name="124" class="C124" element="C" mass="12.01078"/>
  <Type name="125" class="H125" element="H" mass="1.007947"/>
  <Type name="126" class="C126" element="C" mass="12.01078"/>
  <Type name="127" class="O127" element="O" mass="15.99943"/>
  <Type name="128" class="O128" element="O" mass="15.99943"/>
  <Type name="129" class="C129" element="C" mass="12.01078"/>
  <Type name="130" class="O130" element="O" mass="15.99943"/>
  <Type name="131" class="N131" element="N" mass="14.00672"/>
  <Type name="132" class="H132" element="H" mass="1.007947"/>
  <Type name="133" class="C133" element="C" mass="12.01078"/>
  <Type name="134" class="H134" element="H" mass="1.007947"/>
  <Type name="135" class="C135" element="C" mass="12.01078"/>
  <Type name="136" class="O136" element="O" mass="15.99943"/>
  <Type name="137" class="N137" element="N" mass="14.00672"/>
  <Type name="138" class="H138" element="H" mass="1.007947"/>
  <Type name="139" class="C139" element="C" mass="12.01078"/>
  <Type name="140" class="H140" element="H" mass="1.007947"/>
  <Type name="141" class="C141" element="C" mass="12.01078"/>
  <Type name="142" class="H142" element="H" mass="1.007947"/>
  <Type name="143" class="C143" element="C" mass="12.01078"/>
  <Type name="144" class="N144" element="N" mass="14.00672"/>
  <Type name="145" class="H145" element="H" mass="1.007947"/>
  <Type name="146" class="C146" element="C" mass="12.01078"/>
  <Type name="147" class="H147" element="H" mass="1.007947"/>
  <Type name="148" class="N148" element="N" mass="14.00672"/>
  <Type name="149" class="C149" element="C" mass="12.01078"/>
  <Type name="150" class="H150" element="H" mass="1.007947"/>
  <Type name="151" class="C151" element="C" mass="12.01078"/>
  <Type name="152" class="O152" element="O" mass="15.99943"/>
  <Type name="153" class="N153" element="N" mass="14.00672"/>
  <Type name="154" class="H154" element="H" mass="1.007947"/>
  <Type name="155" class="C155" element="C" mass="12.01078"/>
  <Type name="156" class="H156" element="H" mass="1.007947"/>
  <Type name="157" class="C157" element="C" mass="12.01078"/>
  <Type name="158" class="H158" element="H" mass="1.007947"/>
  <Type name="159" class="C159" element="C" mass="12.01078"/>
  <Type name="160" class="N160" element="N" mass="14.00672"/>
  <Type name="161" class="C161" element="C" mass="12.01078"/>
  <Type name="162" class="H162" element="H" mass="1.007947"/>
  <Type name="163" class="N163" element="N" mass="14.00672"/>
  <Type name="164" class="H164" element="H" mass="1.007947"/>
  <Type name="165" class="C165" element="C" mass="12.01078"/>
  <Type name="166" class="H166" element="H" mass="1.007947"/>
  <Type name="167" class="C167" element="C" mass="12.01078"/>
  <Type name="168" class="O168" element="O" mass="15.99943"/>
  <Type name="169" class="N169" element="N" mass="14.00672"/>
  <Type name="170" class="H170" element="H" mass="1.007947"/>
  <Type name="171" class="C171" element="C" mass="12.01078"/>
  <Type name="172" class="H172" element="H" mass="1.007947"/>
  <Type name="173" class="C173" element="C" mass="12.01078"/>
  <Type name="174" class="H174" element="H" mass="1.007947"/>
  <Type name="175" class="C175" element="C" mass="12.01078"/>
  <Type name="176" class="N176" element="N" mass="14.00672"/>
  <Type name="177" class="H177" element="H" mass="1.007947"/>
  <Type name="178" class="C178" element="C" mass="12.01078"/>
  <Type name="179" class="H179" element="H" mass="1.007947"/>
  <Type name="180" class="N180" element="N" mass="14.00672"/>
  <Type name="181" class="H181" element="H" mass="1.007947"/>
  <Type name="182" class="C182" element="C" mass="12.01078"/>
  <Type name="183" class="H183" element="H" mass="1.007947"/>
  <Type name="184" class="C184" element="C" mass="12.01078"/>
  <Type name="185" class="O185" element="O" mass="15.99943"/>
  <Type name="186" class="N186" element="N" mass="14.00672"/>
  <Type name="187" class="H187" element="H" mass="1.007947"/>
  <Type name="188" class="C188" element="C" mass="12.01078"/>
  <Type name="189" class="H189" element="H" mass="1.007947"/>
  <Type name="190" class="C190" element="C" mass="12.01078"/>
  <Type name="191" class="H191" element="H" mass="1.007947"/>
  <Type name="192" class="C192" element="C" mass="12.01078"/>
  <Type name="193" class="H193" element="H" mass="1.007947"/>
  <Type name="194" class="C194" element="C" mass="12.01078"/>
  <Type name="195" class="H195" element="H" mass="1.007947"/>
  <Type name="196" class="C196" element="C" mass="12.01078"/>
  <Type name="197" class="H197" element="H" mass="1.007947"/>
  <Type name="198" class="C198" element="C" mass="12.01078"/>
  <Type name="199" class="O199" element="O" mass="15.99943"/>
  <Type name="200" class="N200" element="N" mass="14.00672"/>
  <Type name="201" class="H201" element="H" mass="1.007947"/>
  <Type name="202" class="C202" element="C" mass="12.01078"/>
  <Type name="203" class="H203" element="H" mass="1.007947"/>
  <Type name="204" class="C204" element="C" mass="12.01078"/>
  <Type name="205" class="H205" element="H" mass="1.007947"/>
  <Type name="206" class="C206" element="C" mass="12.01078"/>
  <Type name="207" class="H207" element="H" mass="1.007947"/>
  <Type name="208" class="C208" element="C" mass="12.01078"/>
  <Type name="209" class="H209" element="H" mass="1.007947"/>
  <Type name="210" class="C210" element="C" mass="12.01078"/>
  <Type name="211" class="H211" element="H" mass="1.007947"/>
  <Type name="212" class="C212" element="C" mass="12.01078"/>
  <Type name="213" class="O213" element="O" mass="15.99943"/>
  <Type name="214" class="N214" element="N" mass="14.00672"/>
  <Type name="215" class="H215" element="H" mass="1.007947"/>
  <Type name="216" class="C216" element="C" mass="12.01078"/>
  <Type name="217" class="H217" element="H" mass="1.007947"/>
  <Type name="218" class="C218" element="C" mass="12.01078"/>
  <Type name="219" class="H219" element="H" mass="1.007947"/>
  <Type name="220" class="C220" element="C" mass="12.01078"/>
  <Type name="221" class="H221" element="H" mass="1.007947"/>
  <Type name="222" class="C222" element="C" mass="12.01078"/>
  <Type name="223" class="H223" element="H" mass="1.007947"/>
  <Type name="224" class="C224" element="C" mass="12.01078"/>
  <Type name="225" class="H225" element="H" mass="1.007947"/>
  <Type name="226" class="N226" element="N" mass="14.00672"/>
  <Type name="227" class="H227" element="H" mass="1.007947"/>
  <Type name="228" class="C228" element="C" mass="12.01078"/>
  <Type name="229" class="O229" element="O" mass="15.99943"/>
  <Type name="230" class="N230" element="N" mass="14.00672"/>
  <Type name="231" class="H231" element="H" mass="1.007947"/>
  <Type name="232" class="C232" element="C" mass="12.01078"/>
  <Type name="233" class="H233" element="H" mass="1.007947"/>
  <Type name="234" class="C234" element="C" mass="12.01078"/>
  <Type name="235" class="H235" element="H" mass="1.007947"/>
  <Type name="236" class="C236" element="C" mass="12.01078"/>
  <Type name="237" class="H237" element="H" mass="1.007947"/>
  <Type name="238" class="C238" element="C" mass="12.01078"/>
  <Type name="239" class="H239" element="H" mass="1.007947"/>
  <Type name="240" class="C240" element="C" mass="12.01078"/>
  <Type name="241" class="H241" element="H" mass="1.007947"/>
  <Type name="242" class="N242" element="N" mass="14.00672"/>
  <Type name="243" class="H243" element="H" mass="1.007947"/>
  <Type name="244" class="C244" element="C" mass="12.01078"/>
  <Type name="245" class="O245" element="O" mass="15.99943"/>
  <Type name="246" class="N246" element="N" mass="14.00672"/>
  <Type name="247" class="H247" element="H" mass="1.007947"/>
  <Type name="248" class="C248" element="C" mass="12.01078"/>
  <Type name="249" class="H249" element="H" mass="1.007947"/>
  <Type name="250" class="C250" element="C" mass="12.01078"/>
  <Type name="251" class="H251" element="H" mass="1.007947"/>
  <Type name="252" class="C252" element="C" mass="12.01078"/>
  <Type name="253" class="H253" element="H" mass="1.007947"/>
  <Type name="254" class="S254" element="S" mass="32.0655"/>
  <Type name="255" class="C255" element="C" mass="12.01078"/>
  <Type name="256" class="H256" element="H" mass="1.007947"/>
  <Type name="257" class="C257" element="C" mass="12.01078"/>
  <Type name="258" class="O258" element="O" mass="15.99943"/>
  <Type name="259" class="N259" element="N" mass="14.00672"/>
  <Type name="260" class="H260" element="H" mass="1.007947"/>
  <Type name="261" class="C261" element="C" mass="12.01078"/>
  <Type name="262" class="H262" element="H" mass="1.007947"/>
  <Type name="263" class="C263" element="C" mass="12.01078"/>
  <Type name="264" class="H264" element="H" mass="1.007947"/>
  <Type name="265" class="C265" element="C" mass="12.01078"/>
  <Type name="266" class="C266" element="C" mass="12.01078"/>
  <Type name="267" class="H267" element="H" mass="1.007947"/>
  <Type name="268" class="C268" element="C" mass="12.01078"/>
  <Type name="269" class="H269" element="H" mass="1.007947"/>
  <Type name="270" class="C270" element="C" mass="12.01078"/>
  <Type name="271" class="H271" element="H" mass="1.007947"/>
  <Type name="272" class="C272" element="C" mass="12.01078"/>
  <Type name="273" class="H273" element="H" mass="1.007947"/>
  <Type name="274" class="C274" element="C" mass="12.01078"/>
  <Type name="275" class="H275" element="H" mass="1.007947"/>
  <Type name="276" class="C276" element="C" mass="12.01078"/>
  <Type name="277" class="O277" element="O" mass="15.99943"/>
  <Type name="278" class="N278" element="N" mass="14.00672"/>
  <Type name="279" class="C279" element="C" mass="12.01078"/>
  <Type name="280" class="H280" element="H" mass="1.007947"/>
  <Type name="281" class="C281" element="C" mass="12.01078"/>
  <Type name="282" class="H282" element="H" mass="1.007947"/>
  <Type name="283" class="C283" element="C" mass="12.01078"/>
  <Type name="284" class="H284" element="H" mass="1.007947"/>
  <Type name="285" class="C285" element="C" mass="12.01078"/>
  <Type name="286" class="H286" element="H" mass="1.007947"/>
  <Type name="287" class="C287" element="C" mass="12.01078"/>
  <Type name="288" class="O288" element="O" mass="15.99943"/>
  <Type name="289" class="N289" element="N" mass="14.00672"/>
  <Type name="290" class="H290" element="H" mass="1.007947"/>
  <Type name="291" class="C291" element="C" mass="12.01078"/>
  <Type name="292" class="H292" element="H" mass="1.007947"/>
  <Type name="293" class="C293" element="C" mass="12.01078"/>
  <Type name="294" class="H294" element="H" mass="1.007947"/>
  <Type name="295" class="O295" element="O" mass="15.99943"/>
  <Type name="296" class="H296" element="H" mass="1.007947"/>
  <Type name="297" class="C297" element="C" mass="12.01078"/>
  <Type name="298" class="O298" element="O" mass="15.99943"/>
  <Type name="299" class="N299" element="N" mass="14.00672"/>
  <Type name="300" class="H300" element="H" mass="1.007947"/>
  <Type name="301" class="C301" element="C" mass="12.01078"/>
  <Type name="302" class="H302" element="H" mass="1.007947"/>
  <Type name="303" class="C303" element="C" mass="12.01078"/>
  <Type name="304" class="H304" element="H" mass="1.007947"/>
  <Type name="305" class="C305" element="C" mass="12.01078"/>
  <Type name="306" class="H306" element="H" mass="1.007947"/>
  <Type name="307" class="O307" element="O" mass="15.99943"/>
  <Type name="308" class="H308" element="H" mass="1.007947"/>
  <Type name="309" class="C309" element="C" mass="12.01078"/>
  <Type name="310" class="O310" element="O" mass="15.99943"/>
  <Type name="311" class="N311" element="N" mass="14.00672"/>
  <Type name="312" class="H312" element="H" mass="1.007947"/>
  <Type name="313" class="C313" element="C" mass="12.01078"/>
  <Type name="314" class="H314" element="H" mass="1.007947"/>
  <Type name="315" class="C315" element="C" mass="12.01078"/>
  <Type name="316" class="H316" element="H" mass="1.007947"/>
  <Type name="317" class="C317" element="C" mass="12.01078"/>
  <Type name="318" class="C318" element="C" mass="12.01078"/>
  <Type name="319" class="H319" element="H" mass="1.007947"/>
  <Type name="320" class="N320" element="N" mass="14.00672"/>
  <Type name="321" class="H321" element="H" mass="1.007947"/>
  <Type name="322" class="C322" element="C" mass="12.01078"/>
  <Type name="323" class="C323" element="C" mass="12.01078"/>
  <Type name="324" class="H324" element="H" mass="1.007947"/>
  <Type name="325" class="C325" element="C" mass="12.01078"/>
  <Type name="326" class="H326" element="H" mass="1.007947"/>
  <Type name="327" class="C327" element="C" mass="12.01078"/>
  <Type name="328" class="H328" element="H" mass="1.007947"/>
  <Type name="329" class="C329" element="C" mass="12.01078"/>
  <Type name="330" class="H330" element="H" mass="1.007947"/>
  <Type name="331" class="C331" element="C" mass="12.01078"/>
  <Type name="332" class="C332" element="C" mass="12.01078"/>
  <Type name="333" class="O333" element="O" mass="15.99943"/>
  <Type name="334" class="N334" element="N" mass="14.00672"/>
  <Type name="335" class="H335" element="H" mass="1.007947"/>
  <Type name="336" class="C336" element="C" mass="12.01078"/>
  <Type name="337" class="H337" element="H" mass="1.007947"/>
  <Type name="338" class="C338" element="C" mass="12.01078"/>
  <Type name="339" class="H339" element="H" mass="1.007947"/>
  <Type name="340" class="C340" element="C" mass="12.01078"/>
  <Type name="341" class="C341" element="C" mass="12.01078"/>
  <Type name="342" class="H342" element="H" mass="1.007947"/>
  <Type name="343" class="C343" element="C" mass="12.01078"/>
  <Type name="344" class="H344" element="H" mass="1.007947"/>
  <Type name="345" class="C345" element="C" mass="12.01078"/>
  <Type name="346" class="O346" element="O" mass="15.99943"/>
  <Type name="347" class="H347" element="H" mass="1.007947"/>
  <Type name="348" class="C348" element="C" mass="12.01078"/>
  <Type name="349" class="H349" element="H" mass="1.007947"/>
  <Type name="350" class="C350" element="C" mass="12.01078"/>
  <Type name="351" class="H351" element="H" mass="1.007947"/>
  <Type name="352" class="C352" element="C" mass="12.01078"/>
  <Type name="353" class="O353" element="O" mass="15.99943"/>
  <Type name="354" class="N354" element="N" mass="14.00672"/>
  <Type name="355" class="H355" element="H" mass="1.007947"/>
  <Type name="356" class="C356" element="C" mass="12.01078"/>
  <Type name="357" class="H357" element="H" mass="1.007947"/>
  <Type name="358" class="C358" element="C" mass="12.01078"/>
  <Type name="359" class="H359" element="H" mass="1.007947"/>
  <Type name="360" class="C360" element="C" mass="12.01078"/>
  <Type name="361" class="H361" element="H" mass="1.007947"/>
  <Type name="362" class="C362" element="C" mass="12.01078"/>
  <Type name="363" class="H363" element="H" mass="1.007947"/>
  <Type name="364" class="C364" element="C" mass="12.01078"/>
  <Type name="365" class="O365" element="O" mass="15.99943"/>
  <Type name="366" class="N366" element="N" mass="14.00672"/>
  <Type name="367" class="H367" element="H" mass="1.007947"/>
  <Type name="368" class="C368" element="C" mass="12.01078"/>
  <Type name="369" class="H369" element="H" mass="1.007947"/>
  <Type name="370" class="C370" element="C" mass="12.01078"/>
  <Type name="371" class="H371" element="H" mass="1.007947"/>
  <Type name="372" class="C372" element="C" mass="12.01078"/>
  <Type name="373" class="O373" element="O" mass="15.99943"/>
  <Type name="374" class="O374" element="O" mass="15.99943"/>
  <Type name="375" class="N375" element="N" mass="14.00672"/>
  <Type name="376" class="H376" element="H" mass="1.007947"/>
  <Type name="377" class="C377" element="C" mass="12.01078"/>
  <Type name="378" class="H378" element="H" mass="1.007947"/>
  <Type name="379" class="C379" element="C" mass="12.01078"/>
  <Type name="380" class="H380" element="H" mass="1.007947"/>
  <Type name="381" class="C381" element="C" mass="12.01078"/>
  <Type name="382" class="H382" element="H" mass="1.007947"/>
  <Type name="383" class="C383" element="C" mass="12.01078"/>
  <Type name="384" class="H384" element="H" mass="1.007947"/>
  <Type name="385" class="N385" element="N" mass="14.00672"/>
  <Type name="386" class="H386" element="H" mass="1.007947"/>
  <Type name="387" class="C387" element="C" mass="12.01078"/>
  <Type name="388" class="N388" element="N" mass="14.00672"/>
  <Type name="389" class="H389" element="H" mass="1.007947"/>
  <Type name="390" class="N390" element="N" mass="14.00672"/>
  <Type name="391" class="H391" element="H" mass="1.007947"/>
  <Type name="392" class="C392" element="C" mass="12.01078"/>
  <Type name="393" class="O393" element="O" mass="15.99943"/>
  <Type name="394" class="O394" element="O" mass="15.99943"/>
  <Type name="395" class="N395" element="N" mass="14.00672"/>
  <Type name="396" class="H396" element="H" mass="1.007947"/>
  <Type name="397" class="C397" element="C" mass="12.01078"/>
  <Type name="398" class="H398" element="H" mass="1.007947"/>
  <Type name="399" class="C399" element="C" mass="12.01078"/>
  <Type name="400" class="H400" element="H" mass="1.007947"/>
  <Type name="401" class="C401" element="C" mass="12.01078"/>
  <Type name="402" class="O402" element="O" mass="15.99943"/>
  <Type name="403" class="N403" element="N" mass="14.00672"/>
  <Type name="404" class="H404" element="H" mass="1.007947"/>
  <Type name="405" class="C405" element="C" mass="12.01078"/>
  <Type name="406" class="O406" element="O" mass="15.99943"/>
  <Type name="407" class="O407" element="O" mass="15.99943"/>
  <Type name="408" class="N408" element="N" mass="14.00672"/>
  <Type name="409" class="H409" element="H" mass="1.007947"/>
  <Type name="410" class="C410" element="C" mass="12.01078"/>
  <Type name="411" class="H411" element="H" mass="1.007947"/>
  <Type name="412" class="C412" element="C" mass="12.01078"/>
  <Type name="413" class="H413" element="H" mass="1.007947"/>
  <Type name="414" class="C414" element="C" mass="12.01078"/>
  <Type name="415" class="O415" element="O" mass="15.99943"/>
  <Type name="416" class="O416" element="O" mass="15.99943"/>
  <Type name="417" class="C417" element="C" mass="12.01078"/>
  <Type name="418" class="O418" element="O" mass="15.99943"/>
  <Type name="419" class="O419" element="O" mass="15.99943"/>
  <Type name="420" class="N420" element="N" mass="14.00672"/>
  <Type name="421" class="H421" element="H" mass="1.007947"/>
  <Type name="422" class="C422" element="C" mass="12.01078"/>
  <Type name="423" class="H423" element="H" mass="1.007947"/>
  <Type name="424" class="C424" element="C" mass="12.01078"/>
  <Type name="425" class="H425" element="H" mass="1.007947"/>
  <Type name="426" class="S426" element="S" mass="32.0655"/>
  <Type name="427" class="H427" element="H" mass="1.007947"/>
  <Type name="428" class="C428" element="C" mass="12.01078"/>
  <Type name="429" class="O429" element="O" mass="15.99943"/>
  <Type name="430" class="O430" element="O" mass="15.99943"/>
  <Type name="431" class="N431" element="N" mass="14.00672"/>
  <Type name="432" class="H432" element="H" mass="1.007947"/>
  <Type name="433" class="C433" element="C" mass="12.01078"/>
  <Type name="434" class="H434" element="H" mass="1.007947"/>
  <Type name="435" class="C435" element="C" mass="12.01078"/>
  <Type name="436" class="H436" element="H" mass="1.007947"/>
  <Type name="437" class="S437" element="S" mass="32.0655"/>
  <Type name="438" class="C438" element="C" mass="12.01078"/>
  <Type name="439" class="O439" element="O" mass="15.99943"/>
  <Type name="440" class="O440" element="O" mass="15.99943"/>
  <Type name="441" class="N441" element="N" mass="14.00672"/>
  <Type name="442" class="H442" element="H" mass="1.007947"/>
  <Type name="443" class="C443" element="C" mass="12.01078"/>
  <Type name="444" class="H444" element="H" mass="1.007947"/>
  <Type name="445" class="C445" element="C" mass="12.01078"/>
  <Type name="446" class="H446" element="H" mass="1.007947"/>
  <Type name="447" class="C447" element="C" mass="12.01078"/>
  <Type name="448" class="H448" element="H" mass="1.007947"/>
  <Type name="449" class="C449" element="C" mass="12.01078"/>
  <Type name="450" class="O450" element="O" mass="15.99943"/>
  <Type name="451" class="N451" element="N" mass="14.00672"/>
  <Type name="452" class="H452" element="H" mass="1.007947"/>
  <Type name="453" class="C453" element="C" mass="12.01078"/>
  <Type name="454" class="O454" element="O" mass="15.99943"/>
  <Type name="455" class="O455" element="O" mass="15.99943"/>
  <Type name="456" class="N456" element="N" mass="14.00672"/>
  <Type name="457" class="H457" element="H" mass="1.007947"/>
  <Type name="458" class="C458" element="C" mass="12.01078"/>
  <Type name="459" class="H459" element="H" mass="1.007947"/>
  <Type name="460" class="C460" element="C" mass="12.01078"/>
  <Type name="461" class="H461" element="H" mass="1.007947"/>
  <Type name="462" class="C462" element="C" mass="12.01078"/>
  <Type name="463" class="H463" element="H" mass="1.007947"/>
  <Type name="464" class="C464" element="C" mass="12.01078"/>
  <Type name="465" class="O465" element="O" mass="15.99943"/>
  <Type name="466" class="O466" element="O" mass="15.99943"/>
  <Type name="467" class="C467" element="C" mass="12.01078"/>
  <Type name="468" class="O468" element="O" mass="15.99943"/>
  <Type name="469" class="O469" element="O" mass="15.99943"/>
  <Type name="470" class="N470" element="N" mass="14.00672"/>
  <Type name="471" class="H471" element="H" mass="1.007947"/>
  <Type name="472" class="C472" element="C" mass="12.01078"/>
  <Type name="473" class="H473" element="H" mass="1.007947"/>
  <Type name="474" class="C474" element="C" mass="12.01078"/>
  <Type name="475" class="O475" element="O" mass="15.99943"/>
  <Type name="476" class="O476" element="O" mass="15.99943"/>
  <Type name="477" class="N477" element="N" mass="14.00672"/>
  <Type name="478" class="H478" element="H" mass="1.007947"/>
  <Type name="479" class="C479" element="C" mass="12.01078"/>
  <Type name="480" class="H480" element="H" mass="1.007947"/>
  <Type name="481" class="C481" element="C" mass="12.01078"/>
  <Type name="482" class="H482" element="H" mass="1.007947"/>
  <Type name="483" class="C483" element="C" mass="12.01078"/>
  <Type name="484" class="N484" element="N" mass="14.00672"/>
  <Type name="485" class="H485" element="H" mass="1.007947"/>
  <Type name="486" class="C486" element="C" mass="12.01078"/>
  <Type name="487" class="H487" element="H" mass="1.007947"/>
  <Type name="488" class="N488" element="N" mass="14.00672"/>
  <Type name="489" class="C489" element="C" mass="12.01078"/>
  <Type name="490" class="H490" element="H" mass="1.007947"/>
  <Type name="491" class="C491" element="C" mass="12.01078"/>
  <Type name="492" class="O492" element="O" mass="15.99943"/>
  <Type name="493" class="O493" element="O" mass="15.99943"/>
  <Type name="494" class="N494" element="N" mass="14.00672"/>
  <Type name="495" class="H495" element="H" mass="1.007947"/>
  <Type name="496" class="C496" element="C" mass="12.01078"/>
  <Type name="497" class="H497" element="H" mass="1.007947"/>
  <Type name="498" class="C498" element="C" mass="12.01078"/>
  <Type name="499" class="H499" element="H" mass="1.007947"/>
  <Type name="500" class="C500" element="C" mass="12.01078"/>
  <Type name="501" class="N501" element="N" mass="14.00672"/>
  <Type name="502" class="C502" element="C" mass="12.01078"/>
  <Type name="503" class="H503" element="H" mass="1.007947"/>
  <Type name="504" class="N504" element="N" mass="14.00672"/>
  <Type name="505" class="H505" element="H" mass="1.007947"/>
  <Type name="506" class="C506" element="C" mass="12.01078"/>
  <Type name="507" class="H507" element="H" mass="1.007947"/>
  <Type name="508" class="C508" element="C" mass="12.01078"/>
  <Type name="509" class="O509" element="O" mass="15.99943"/>
  <Type name="510" class="O510" element="O" mass="15.99943"/>
  <Type name="511" class="N511" element="N" mass="14.00672"/>
  <Type name="512" class="H512" element="H" mass="1.007947"/>
  <Type name="513" class="C513" element="C" mass="12.01078"/>
  <Type name="514" class="H514" element="H" mass="1.007947"/>
  <Type name="515" class="C515" element="C" mass="12.01078"/>
  <Type name="516" class="H516" element="H" mass="1.007947"/>
  <Type name="517" class="C517" element="C" mass="12.01078"/>
  <Type name="518" class="N518" element="N" mass="14.00672"/>
  <Type name="519" class="H519" element="H" mass="1.007947"/>
  <Type name="520" class="C520" element="C" mass="12.01078"/>
  <Type name="521" class="H521" element="H" mass="1.007947"/>
  <Type name="522" class="N522" element="N" mass="14.00672"/>
  <Type name="523" class="H523" element="H" mass="1.007947"/>
  <Type name="524" class="C524" element="C" mass="12.01078"/>
  <Type name="525" class="H525" element="H" mass="1.007947"/>
  <Type name="526" class="C" element="C" mass="12.01078"/>
  <Type name="527" class="O2" element="O" mass="15.99943"/>
  <Type name="528" class="O2" element="O" mass="15.99943"/>
  <Type name="529" class="N" element="N" mass="14.00672"/>
  <Type name="530" class="H" element="H" mass="1.007947"/>
  <Type name="531" class="CT" element="C" mass="12.01078"/>
  <Type name="532" class="H1" element="H" mass="1.007947"/>
  <Type name="533" class="CT" element="C" mass="12.01078"/>
  <Type name="534" class="HC" element="H" mass="1.007947"/>
  <Type name="535" class="CT" element="C" mass="12.01078"/>
  <Type name="536" class="HC" element="H" mass="1.007947"/>
  <Type name="537" class="CT" element="C" mass="12.01078"/>
  <Type name="538" class="HC" element="H" mass="1.007947"/>
  <Type name="539" class="CT" element="C" mass="12.01078"/>
  <Type name="540" class="HC" element="H" mass="1.007947"/>
  <Type name="541" class="C" element="C" mass="12.01078"/>
  <Type name="542" class="O2" element="O" mass="15.99943"/>
  <Type name="543" class="O2" element="O" mass="15.99943"/>
  <Type name="544" class="N" element="N" mass="14.00672"/>
  <Type name="545" class="H" element="H" mass="1.007947"/>
  <Type name="546" class="CT" element="C" mass="12.01078"/>
  <Type name="547" class="H1" element="H" mass="1.007947"/>
  <Type name="548" class="CT" element="C" mass="12.01078"/>
  <Type name="549" class="HC" element="H" mass="1.007947"/>
  <Type name="550" class="CT" element="C" mass="12.01078"/>
  <Type name="551" class="HC" element="H" mass="1.007947"/>
  <Type name="552" class="CT" element="C" mass="12.01078"/>
  <Type name="553" class="HC" element="H" mass="1.007947"/>
  <Type name="554" class="CT" element="C" mass="12.01078"/>
  <Type name="555" class="HC" element="H" mass="1.007947"/>
  <Type name="556" class="C" element="C" mass="12.01078"/>
  <Type name="557" class="O2" element="O" mass="15.99943"/>
  <Type name="558" class="O2" element="O" mass="15.99943"/>
  <Type name="559" class="N" element="N" mass="14.00672"/>
  <Type name="560" class="H" element="H" mass="1.007947"/>
  <Type name="561" class="CT" element="C" mass="12.01078"/>
  <Type name="562" class="H1" element="H" mass="1.007947"/>
  <Type name="563" class="CT" element="C" mass="12.01078"/>
  <Type name="564" class="HC" element="H" mass="1.007947"/>
  <Type name="565" class="CT" element="C" mass="12.01078"/>
  <Type name="566" class="HC" element="H" mass="1.007947"/>
  <Type name="567" class="CT" element="C" mass="12.01078"/>
  <Type name="568" class="HC" element="H" mass="1.007947"/>
  <Type name="569" class="CT" element="C" mass="12.01078"/>
  <Type name="570" class="HP" element="H" mass="1.007947"/>
  <Type name="571" class="N3" element="N" mass="14.00672"/>
  <Type name="572" class="H" element="H" mass="1.007947"/>
  <Type name="573" class="C" element="C" mass="12.01078"/>
  <Type name="574" class="O2" element="O" mass="15.99943"/>
  <Type name="575" class="O2" element="O" mass="15.99943"/>
  <Type name="576" class="N" element="N" mass="14.00672"/>
  <Type name="577" class="H" element="H" mass="1.007947"/>
  <Type name="578" class="CT" element="C" mass="12.01078"/>
  <Type name="579" class="H1" element="H" mass="1.007947"/>
  <Type name="580" class="CT" element="C" mass="12.01078"/>
  <Type name="581" class="HC" element="H" mass="1.007947"/>
  <Type name="582" class="CT" element="C" mass="12.01078"/>
  <Type name="583" class="H1" element="H" mass="1.007947"/>
  <Type name="584" class="S" element="S" mass="32.0655"/>
  <Type name="585" class="CT" element="C" mass="12.01078"/>
  <Type name="586" class="H1" element="H" mass="1.007947"/>
  <Type name="587" class="C" element="C" mass="12.01078"/>
  <Type name="588" class="O2" element="O" mass="15.99943"/>
  <Type name="589" class="O2" element="O" mass="15.99943"/>
  <Type name="590" class="N" element="N" mass="14.00672"/>
  <Type name="591" class="H" element="H" mass="1.007947"/>
  <Type name="592" class="CT" element="C" mass="12.01078"/>
  <Type name="593" class="H1" element="H" mass="1.007947"/>
  <Type name="594" class="CT" element="C" mass="12.01078"/>
  <Type name="595" class="HC" element="H" mass="1.007947"/>
  <Type name="596" class="CA" element="C" mass="12.01078"/>
  <Type name="597" class="CA" element="C" mass="12.01078"/>
  <Type name="598" class="HA" element="H" mass="1.007947"/>
  <Type name="599" class="CA" element="C" mass="12.01078"/>
  <Type name="600" class="HA" element="H" mass="1.007947"/>
  <Type name="601" class="CA" element="C" mass="12.01078"/>
  <Type name="602" class="HA" element="H" mass="1.007947"/>
  <Type name="603" class="CA" element="C" mass="12.01078"/>
  <Type name="604" class="HA" element="H" mass="1.007947"/>
  <Type name="605" class="CA" element="C" mass="12.01078"/>
  <Type name="606" class="HA" element="H" mass="1.007947"/>
  <Type name="607" class="C" element="C" mass="12.01078"/>
  <Type name="608" class="O2" element="O" mass="15.99943"/>
  <Type name="609" class="O2" element="O" mass="15.99943"/>
  <Type name="610" class="N" element="N" mass="14.00672"/>
  <Type name="611" class="CT" element="C" mass="12.01078"/>
  <Type name="612" class="H1" element="H" mass="1.007947"/>
  <Type name="613" class="CT" element="C" mass="12.01078"/>
  <Type name="614" class="HC" element="H" mass="1.007947"/>
  <Type name="615" class="CT" element="C" mass="12.01078"/>
  <Type name="616" class="HC" element="H" mass="1.007947"/>
  <Type name="617" class="CT" element="C" mass="12.01078"/>
  <Type name="618" class="H1" element="H" mass="1.007947"/>
  <Type name="619" class="C" element="C" mass="12.01078"/>
  <Type name="620" class="O2" element="O" mass="15.99943"/>
  <Type name="621" class="O2" element="O" mass="15.99943"/>
  <Type name="622" class="N" element="N" mass="14.00672"/>
  <Type name="623" class="H" element="H" mass="1.007947"/>
  <Type name="624" class="CT" element="C" mass="12.01078"/>
  <Type name="625" class="H1" element="H" mass="1.007947"/>
  <Type name="626" class="CT" element="C" mass="12.01078"/>
  <Type name="627" class="H1" element="H" mass="1.007947"/>
  <Type name="628" class="OH" element="O" mass="15.99943"/>
  <Type name="629" class="HO" element="H" mass="1.007947"/>
  <Type name="630" class="C" element="C" mass="12.01078"/>
  <Type name="631" class="O2" element="O" mass="15.99943"/>
  <Type name="632" class="O2" element="O" mass="15.99943"/>
  <Type name="633" class="N" element="N" mass="14.00672"/>
  <Type name="634" class="H" element="H" mass="1.007947"/>
  <Type name="635" class="CT" element="C" mass="12.01078"/>
  <Type name="636" class="H1" element="H" mass="1.007947"/>
  <Type name="637" class="CT" element="C" mass="12.01078"/>
  <Type name="638" class="H1" element="H" mass="1.007947"/>
  <Type name="639" class="CT" element="C" mass="12.01078"/>
  <Type name="640" class="HC" element="H" mass="1.007947"/>
  <Type name="641" class="OH" element="O" mass="15.99943"/>
  <Type name="642" class="HO" element="H" mass="1.007947"/>
  <Type name="643" class="C" element="C" mass="12.01078"/>
  <Type name="644" class="O2" element="O" mass="15.99943"/>
  <Type name="645" class="O2" element="O" mass="15.99943"/>
  <Type name="646" class="N" element="N" mass="14.00672"/>
  <Type name="647" class="H" element="H" mass="1.007947"/>
  <Type name="648" class="CT" element="C" mass="12.01078"/>
  <Type name="649" class="H1" element="H" mass="1.007947"/>
  <Type name="650" class="CT" element="C" mass="12.01078"/>
  <Type name="651" class="HC" element="H" mass="1.007947"/>
  <Type name="652" class="C*" element="C" mass="12.01078"/>
  <Type name="653" class="CW" element="C" mass="12.01078"/>
  <Type name="654" class="H4" element="H" mass="1.007947"/>
  <Type name="655" class="NA" element="N" mass="14.00672"/>
  <Type name="656" class="H" element="H" mass="1.007947"/>
  <Type name="657" class="CN" element="C" mass="12.01078"/>
  <Type name="658" class="CA" element="C" mass="12.01078"/>
  <Type name="659" class="HA" element="H" mass="1.007947"/>
  <Type name="660" class="CA" element="C" mass="12.01078"/>
  <Type name="661" class="HA" element="H" mass="1.007947"/>
  <Type name="662" class="CA" element="C" mass="12.01078"/>
  <Type name="663" class="HA" element="H" mass="1.007947"/>
  <Type name="664" class="CA" element="C" mass="12.01078"/>
  <Type name="665" class="HA" element="H" mass="1.007947"/>
  <Type name="666" class="CB" element="C" mass="12.01078"/>
  <Type name="667" class="C" element="C" mass="12.01078"/>
  <Type name="668" class="O2" element="O" mass="15.99943"/>
  <Type name="669" class="O2" element="O" mass="15.99943"/>
  <Type name="670" class="N" element="N" mass="14.00672"/>
  <Type name="671" class="H" element="H" mass="1.007947"/>
  <Type name="672" class="CT" element="C" mass="12.01078"/>
  <Type name="673" class="H1" element="H" mass="1.007947"/>
  <Type name="674" class="CT" element="C" mass="12.01078"/>
  <Type name="675" class="HC" element="H" mass="1.007947"/>
  <Type name="676" class="CA" element="C" mass="12.01078"/>
  <Type name="677" class="CA" element="C" mass="12.01078"/>
  <Type name="678" class="HA" element="H" mass="1.007947"/>
  <Type name="679" class="CA" element="C" mass="12.01078"/>
  <Type name="680" class="HA" element="H" mass="1.007947"/>
  <Type name="681" class="C" element="C" mass="12.01078"/>
  <Type name="682" class="OH" element="O" mass="15.99943"/>
  <Type name="683" class="HO" element="H" mass="1.007947"/>
  <Type name="684" class="CA" element="C" mass="12.01078"/>
  <Type name="685" class="HA" element="H" mass="1.007947"/>
  <Type name="686" class="CA" element="C" mass="12.01078"/>
  <Type name="687" class="HA" element="H" mass="1.007947"/>
  <Type name="688" class="C" element="C" mass="12.01078"/>
  <Type name="689" class="O2" element="O" mass="15.99943"/>
  <Type name="690" class="O2" element="O" mass="15.99943"/>
  <Type name="691" class="N" element="N" mass="14.00672"/>
  <Type name="692" class="H" element="H" mass="1.007947"/>
  <Type name="693" class="CT" element="C" mass="12.01078"/>
  <Type name="694" class="H1" element="H" mass="1.007947"/>
  <Type name="695" class="CT" element="C" mass="12.01078"/>
  <Type name="696" class="HC" element="H" mass="1.007947"/>
  <Type name="697" class="CT" element="C" mass="12.01078"/>
  <Type name="698" class="HC" element="H" mass="1.007947"/>
  <Type name="699" class="CT" element="C" mass="12.01078"/>
  <Type name="700" class="HC" element="H" mass="1.007947"/>
  <Type name="701" class="C" element="C" mass="12.01078"/>
  <Type name="702" class="O2" element="O" mass="15.99943"/>
  <Type name="703" class="O2" element="O" mass="15.99943"/>
  <Type name="704" class="N" element="N" mass="14.00672"/>
  <Type name="705" class="H" element="H" mass="1.007947"/>
  <Type name="706" class="N706" element="N" mass="14.00672"/>
  <Type name="707" class="H707" element="H" mass="1.007947"/>
  <Type name="708" class="C708" element="C" mass="12.01078"/>
  <Type name="709" class="H709" element="H" mass="1.007947"/>
  <Type name="710" class="H710" element="H" mass="1.007947"/>
  <Type name="711" class="C711" element="C" mass="12.01078"/>
  <Type name="712" class="C712" element="C" mass="12.01078"/>
  <Type name="713" class="O713" element="O" mass="15.99943"/>
  <Type name="714" class="N714" element="N" mass="14.00672"/>
  <Type name="715" class="H715" element="H" mass="1.007947"/>
  <Type name="716" class="C716" element="C" mass="12.01078"/>
  <Type name="717" class="H717" element="H" mass="1.007947"/>
  <Type name="718" class="C718" element="C" mass="12.01078"/>
  <Type name="719" class="H719" element="H" mass="1.007947"/>
  <Type name="720" class="C720" element="C" mass="12.01078"/>
  <Type name="721" class="O721" element="O" mass="15.99943"/>
  <Type name="722" class="N722" element="N" mass="14.00672"/>
  <Type name="723" class="H723" element="H" mass="1.007947"/>
  <Type name="724" class="C724" element="C" mass="12.01078"/>
  <Type name="725" class="H725" element="H" mass="1.007947"/>
  <Type name="726" class="C726" element="C" mass="12.01078"/>
  <Type name="727" class="HC" element="H" mass="1.007947"/>
  <Type name="728" class="CT" element="C" mass="12.01078"/>
  <Type name="729" class="HC" element="H" mass="1.007947"/>
  <Type name="730" class="CT" element="C" mass="12.01078"/>
  <Type name="731" class="H1" element="H" mass="1.007947"/>
  <Type name="732" class="N2" element="N" mass="14.00672"/>
  <Type name="733" class="H" element="H" mass="1.007947"/>
  <Type name="734" class="CA" element="C" mass="12.01078"/>
  <Type name="735" class="N2" element="N" mass="14.00672"/>
  <Type name="736" class="H" element="H" mass="1.007947"/>
  <Type name="737" class="N2" element="N" mass="14.00672"/>
  <Type name="738" class="H" element="H" mass="1.007947"/>
  <Type name="739" class="C" element="C" mass="12.01078"/>
  <Type name="740" class="O" element="O" mass="15.99943"/>
  <Type name="741" class="N3" element="N" mass="14.00672"/>
  <Type name="742" class="H" element="H" mass="1.007947"/>
  <Type name="743" class="CT" element="C" mass="12.01078"/>
  <Type name="744" class="HP" element="H" mass="1.007947"/>
  <Type name="745" class="CT" element="C" mass="12.01078"/>
  <Type name="746" class="HC" element="H" mass="1.007947"/>
  <Type name="747" class="C" element="C" mass="12.01078"/>
  <Type name="748" class="O" element="O" mass="15.99943"/>
  <Type name="749" class="N" element="N" mass="14.00672"/>
  <Type name="750" class="H" element="H" mass="1.007947"/>
  <Type name="751" class="C" element="C" mass="12.01078"/>
  <Type name="752" class="O" element="O" mass="15.99943"/>
  <Type name="753" class="N3" element="N" mass="14.00672"/>
  <Type name="754" class="H" element="H" mass="1.007947"/>
  <Type name="755" class="CT" element="C" mass="12.01078"/>
  <Type name="756" class="HP" element="H" mass="1.007947"/>
  <Type name="757" class="CT" element="C" mass="12.01078"/>
  <Type name="758" class="HC" element="H" mass="1.007947"/>
  <Type name="759" class="C" element="C" mass="12.01078"/>
  <Type name="760" class="O2" element="O" mass="15.99943"/>
  <Type name="761" class="O2" element="O" mass="15.99943"/>
  <Type name="762" class="C" element="C" mass="12.01078"/>
  <Type name="763" class="O" element="O" mass="15.99943"/>
  <Type name="764" class="N3" element="N" mass="14.00672"/>
  <Type name="765" class="H" element="H" mass="1.007947"/>
  <Type name="766" class="CT" element="C" mass="12.01078"/>
  <Type name="767" class="HP" element="H" mass="1.007947"/>
  <Type name="768" class="CT" element="C" mass="12.01078"/>
  <Type name="769" class="H1" element="H" mass="1.007947"/>
  <Type name="770" class="SH" element="S" mass="32.0655"/>
  <Type name="771" class="HS" element="H" mass="1.007947"/>
  <Type name="772" class="C" element="C" mass="12.01078"/>
  <Type name="773" class="O" element="O" mass="15.99943"/>
  <Type name="774" class="N3" element="N" mass="14.00672"/>
  <Type name="775" class="H" element="H" mass="1.007947"/>
  <Type name="776" class="CT" element="C" mass="12.01078"/>
  <Type name="777" class="HP" element="H" mass="1.007947"/>
  <Type name="778" class="CT" element="C" mass="12.01078"/>
  <Type name="779" class="H1" element="H" mass="1.007947"/>
  <Type name="780" class="S" element="S" mass="32.0655"/>
  <Type name="781" class="C" element="C" mass="12.01078"/>
  <Type name="782" class="O" element="O" mass="15.99943"/>
  <Type name="783" class="N3" element="N" mass="14.00672"/>
  <Type name="784" class="H" element="H" mass="1.007947"/>
  <Type name="785" class="CT" element="C" mass="12.01078"/>
  <Type name="786" class="HP" element="H" mass="1.007947"/>
  <Type name="787" class="CT" element="C" mass="12.01078"/>
  <Type name="788" class="HC" element="H" mass="1.007947"/>
  <Type name="789" class="CT" element="C" mass="12.01078"/>
  <Type name="790" class="HC" element="H" mass="1.007947"/>
  <Type name="791" class="C" element="C" mass="12.01078"/>
  <Type name="792" class="O" element="O" mass="15.99943"/>
  <Type name="793" class="N" element="N" mass="14.00672"/>
  <Type name="794" class="H" element="H" mass="1.007947"/>
  <Type name="795" class="C" element="C" mass="12.01078"/>
  <Type name="796" class="O" element="O" mass="15.99943"/>
  <Type name="797" class="N3" element="N" mass="14.00672"/>
  <Type name="798" class="H" element="H" mass="1.007947"/>
  <Type name="799" class="CT" element="C" mass="12.01078"/>
  <Type name="800" class="HP" element="H" mass="1.007947"/>
  <Type name="801" class="CT" element="C" mass="12.01078"/>
  <Type name="802" class="HC" element="H" mass="1.007947"/>
  <Type name="803" class="CT" element="C" mass="12.01078"/>
  <Type name="804" class="HC" element="H" mass="1.007947"/>
  <Type name="805" class="C" element="C" mass="12.01078"/>
  <Type name="806" class="O2" element="O" mass="15.99943"/>
  <Type name="807" class="O2" element="O" mass="15.99943"/>
  <Type name="808" class="C" element="C" mass="12.01078"/>
  <Type name="809" class="O" element="O" mass="15.99943"/>
  <Type name="810" class="N3" element="N" mass="14.00672"/>
  <Type name="811" class="H" element="H" mass="1.007947"/>
  <Type name="812" class="CT" element="C" mass="12.01078"/>
  <Type name="813" class="HP" element="H" mass="1.007947"/>
  <Type name="814" class="C" element="C" mass="12.01078"/>
  <Type name="815" class="O" element="O" mass="15.99943"/>
  <Type name="816" class="N3" element="N" mass="14.00672"/>
  <Type name="817" class="H" element="H" mass="1.007947"/>
  <Type name="818" class="CT" element="C" mass="12.01078"/>
  <Type name="819" class="HP" element="H" mass="1.007947"/>
  <Type name="820" class="CT" element="C" mass="12.01078"/>
  <Type name="821" class="HC" element="H" mass="1.007947"/>
  <Type name="822" class="CC" element="C" mass="12.01078"/>
  <Type name="823" class="NA" element="N" mass="14.00672"/>
  <Type name="824" class="H" element="H" mass="1.007947"/>
  <Type name="825" class="CR" element="C" mass="12.01078"/>
  <Type name="826" class="H5" element="H" mass="1.007947"/>
  <Type name="827" class="NB" element="N" mass="14.00672"/>
  <Type name="828" class="CV" element="C" mass="12.01078"/>
  <Type name="829" class="H4" element="H" mass="1.007947"/>
  <Type name="830" class="C" element="C" mass="12.01078"/>
  <Type name="831" class="O" element="O" mass="15.99943"/>
  <Type name="832" class="N3" element="N" mass="14.00672"/>
  <Type name="833" class="H" element="H" mass="1.007947"/>
  <Type name="834" class="CT" element="C" mass="12.01078"/>
  <Type name="835" class="HP" element="H" mass="1.007947"/>
  <Type name="836" class="CT" element="C" mass="12.01078"/>
  <Type name="837" class="HC" element="H" mass="1.007947"/>
  <Type name="838" class="CC" element="C" mass="12.01078"/>
  <Type name="839" class="NB" element="N" mass="14.00672"/>
  <Type name="840" class="CR" element="C" mass="12.01078"/>
  <Type name="841" class="H5" element="H" mass="1.007947"/>
  <Type name="842" class="NA" element="N" mass="14.00672"/>
  <Type name="843" class="H" element="H" mass="1.007947"/>
  <Type name="844" class="CW" element="C" mass="12.01078"/>
  <Type name="845" class="H4" element="H" mass="1.007947"/>
  <Type name="846" class="C" element="C" mass="12.01078"/>
  <Type name="847" class="O" element="O" mass="15.99943"/>
  <Type name="848" class="N3" element="N" mass="14.00672"/>
  <Type name="849" class="H" element="H" mass="1.007947"/>
  <Type name="850" class="CT" element="C" mass="12.01078"/>
  <Type name="851" class="HP" element="H" mass="1.007947"/>
  <Type name="852" class="CT" element="C" mass="12.01078"/>
  <Type name="853" class="HC" element="H" mass="1.007947"/>
  <Type name="854" class="CC" element="C" mass="12.01078"/>
  <Type name="855" class="NA" element="N" mass="14.00672"/>
  <Type name="856" class="H" element="H" mass="1.007947"/>
  <Type name="857" class="CR" element="C" mass="12.01078"/>
  <Type name="858" class="H5" element="H" mass="1.007947"/>
  <Type name="859" class="NA" element="N" mass="14.00672"/>
  <Type name="860" class="H" element="H" mass="1.007947"/>
  <Type name="861" class="CW" element="C" mass="12.01078"/>
  <Type name="862" class="H4" element="H" mass="1.007947"/>
  <Type name="863" class="C" element="C" mass="12.01078"/>
  <Type name="864" class="O" element="O" mass="15.99943"/>
  <Type name="865" class="N3" element="N" mass="14.00672"/>
  <Type name="866" class="H" element="H" mass="1.007947"/>
  <Type name="867" class="CT" element="C" mass="12.01078"/>
  <Type name="868" class="HP" element="H" mass="1.007947"/>
  <Type name="869" class="CT" element="C" mass="12.01078"/>
  <Type name="870" class="HC" element="H" mass="1.007947"/>
  <Type name="871" class="CT" element="C" mass="12.01078"/>
  <Type name="872" class="HC" element="H" mass="1.007947"/>
  <Type name="873" class="CT" element="C" mass="12.01078"/>
  <Type name="874" class="HC" element="H" mass="1.007947"/>
  <Type name="875" class="CT" element="C" mass="12.01078"/>
  <Type name="876" class="HC" element="H" mass="1.007947"/>
  <Type name="877" class="C" element="C" mass="12.01078"/>
  <Type name="878" class="O" element="O" mass="15.99943"/>
  <Type name="879" class="N3" element="N" mass="14.00672"/>
  <Type name="880" class="H" element="H" mass="1.007947"/>
  <Type name="881" class="CT" element="C" mass="12.01078"/>
  <Type name="882" class="HP" element="H" mass="1.007947"/>
  <Type name="883" class="CT" element="C" mass="12.01078"/>
  <Type name="884" class="HC" element="H" mass="1.007947"/>
  <Type name="885" class="CT" element="C" mass="12.01078"/>
  <Type name="886" class="HC" element="H" mass="1.007947"/>
  <Type name="887" class="CT" element="C" mass="12.01078"/>
  <Type name="888" class="HC" element="H" mass="1.007947"/>
  <Type name="889" class="CT" element="C" mass="12.01078"/>
  <Type name="890" class="HC" element="H" mass="1.007947"/>
  <Type name="891" class="C" element="C" mass="12.01078"/>
  <Type name="892" class="O" element="O" mass="15.99943"/>
  <Type name="893" class="N3" element="N" mass="14.00672"/>
  <Type name="894" class="H" element="H" mass="1.007947"/>
  <Type name="895" class="CT" element="C" mass="12.01078"/>
  <Type name="896" class="HP" element="H" mass="1.007947"/>
  <Type name="897" class="CT" element="C" mass="12.01078"/>
  <Type name="898" class="HC" element="H" mass="1.007947"/>
  <Type name="899" class="CT" element="C" mass="12.01078"/>
  <Type name="900" class="HC" element="H" mass="1.007947"/>
  <Type name="901" class="CT" element="C" mass="12.01078"/>
  <Type name="902" class="HC" element="H" mass="1.007947"/>
  <Type name="903" class="CT" element="C" mass="12.01078"/>
  <Type name="904" class="HP" element="H" mass="1.007947"/>
  <Type name="905" class="N3" element="N" mass="14.00672"/>
  <Type name="906" class="H" element="H" mass="1.007947"/>
  <Type name="907" class="C" element="C" mass="12.01078"/>
  <Type name="908" class="O" element="O" mass="15.99943"/>
  <Type name="909" class="N3" element="N" mass="14.00672"/>
  <Type name="910" class="H" element="H" mass="1.007947"/>
  <Type name="911" class="CT" element="C" mass="12.01078"/>
  <Type name="912" class="HP" element="H" mass="1.007947"/>
  <Type name="913" class="CT" element="C" mass="12.01078"/>
  <Type name="914" class="HC" element="H" mass="1.007947"/>
  <Type name="915" class="CT" element="C" mass="12.01078"/>
  <Type name="916" class="H1" element="H" mass="1.007947"/>
  <Type name="917" class="S" element="S" mass="32.0655"/>
  <Type name="918" class="CT" element="C" mass="12.01078"/>
  <Type name="919" class="H1" element="H" mass="1.007947"/>
  <Type name="920" class="C" element="C" mass="12.01078"/>
  <Type name="921" class="O" element="O" mass="15.99943"/>
  <Type name="922" class="N3" element="N" mass="14.00672"/>
  <Type name="923" class="H" element="H" mass="1.007947"/>
  <Type name="924" class="CT" element="C" mass="12.01078"/>
  <Type name="925" class="HP" element="H" mass="1.007947"/>
  <Type name="926" class="CT" element="C" mass="12.01078"/>
  <Type name="927" class="HC" element="H" mass="1.007947"/>
  <Type name="928" class="CA" element="C" mass="12.01078"/>
  <Type name="929" class="CA" element="C" mass="12.01078"/>
  <Type name="930" class="HA" element="H" mass="1.007947"/>
  <Type name="931" class="CA" element="C" mass="12.01078"/>
  <Type name="932" class="HA" element="H" mass="1.007947"/>
  <Type name="933" class="CA" element="C" mass="12.01078"/>
  <Type name="934" class="HA" element="H" mass="1.007947"/>
  <Type name="935" class="CA" element="C" mass="12.01078"/>
  <Type name="936" class="HA" element="H" mass="1.007947"/>
  <Type name="937" class="CA" element="C" mass="12.01078"/>
  <Type name="938" class="HA" element="H" mass="1.007947"/>
  <Type name="939" class="C" element="C" mass="12.01078"/>
  <Type name="940" class="O" element="O" mass="15.99943"/>
  <Type name="941" class="N3" element="N" mass="14.00672"/>
  <Type name="942" class="H" element="H" mass="1.007947"/>
  <Type name="943" class="CT" element="C" mass="12.01078"/>
  <Type name="944" class="HP" element="H" mass="1.007947"/>
  <Type name="945" class="CT" element="C" mass="12.01078"/>
  <Type name="946" class="HC" element="H" mass="1.007947"/>
  <Type name="947" class="CT" element="C" mass="12.01078"/>
  <Type name="948" class="HC" element="H" mass="1.007947"/>
  <Type name="949" class="CT" element="C" mass="12.01078"/>
  <Type name="950" class="HP" element="H" mass="1.007947"/>
  <Type name="951" class="C" element="C" mass="12.01078"/>
  <Type name="952" class="O" element="O" mass="15.99943"/>
  <Type name="953" class="N3" element="N" mass="14.00672"/>
  <Type name="954" class="H" element="H" mass="1.007947"/>
  <Type name="955" class="CT" element="C" mass="12.01078"/>
  <Type name="956" class="HP" element="H" mass="1.007947"/>
  <Type name="957" class="CT" element="C" mass="12.01078"/>
  <Type name="958" class="H1" element="H" mass="1.007947"/>
  <Type name="959" class="OH" element="O" mass="15.99943"/>
  <Type name="960" class="HO" element="H" mass="1.007947"/>
  <Type name="961" class="C" element="C" mass="12.01078"/>
  <Type name="962" class="O" element="O" mass="15.99943"/>
  <Type name="963" class="N3" element="N" mass="14.00672"/>
  <Type name="964" class="H" element="H" mass="1.007947"/>
  <Type name="965" class="CT" element="C" mass="12.01078"/>
  <Type name="966" class="HP" element="H" mass="1.007947"/>
  <Type name="967" class="CT" element="C" mass="12.01078"/>
  <Type name="968" class="H1" element="H" mass="1.007947"/>
  <Type name="969" class="CT" element="C" mass="12.01078"/>
  <Type name="970" class="HC" element="H" mass="1.007947"/>
  <Type name="971" class="OH" element="O" mass="15.99943"/>
  <Type name="972" class="HO" element="H" mass="1.007947"/>
  <Type name="973" class="C" element="C" mass="12.01078"/>
  <Type name="974" class="O" element="O" mass="15.99943"/>
  <Type name="975" class="N3" element="N" mass="14.00672"/>
  <Type name="976" class="H" element="H" mass="1.007947"/>
  <Type name="977" class="CT" element="C" mass="12.01078"/>
  <Type name="978" class="HP" element="H" mass="1.007947"/>
  <Type name="979" class="CT" element="C" mass="12.01078"/>
  <Type name="980" class="HC" element="H" mass="1.007947"/>
  <Type name="981" class="C*" element="C" mass="12.01078"/>
  <Type name="982" class="CW" element="C" mass="12.01078"/>
  <Type name="983" class="H4" element="H" mass="1.007947"/>
  <Type name="984" class="NA" element="N" mass="14.00672"/>
  <Type name="985" class="H" element="H" mass="1.007947"/>
  <Type name="986" class="CN" element="C" mass="12.01078"/>
  <Type name="987" class="CA" element="C" mass="12.01078"/>
  <Type name="988" class="HA" element="H" mass="1.007947"/>
  <Type name="989" class="CA" element="C" mass="12.01078"/>
  <Type name="990" class="HA" element="H" mass="1.007947"/>
  <Type name="991" class="CA" element="C" mass="12.01078"/>
  <Type name="992" class="HA" element="H" mass="1.007947"/>
  <Type name="993" class="CA" element="C" mass="12.01078"/>
  <Type name="994" class="HA" element="H" mass="1.007947"/>
  <Type name="995" class="CB" element="C" mass="12.01078"/>
  <Type name="996" class="C" element="C" mass="12.01078"/>
  <Type name="997" class="O" element="O" mass="15.99943"/>
  <Type name="998" class="N3" element="N" mass="14.00672"/>
  <Type name="999" class="H" element="H" mass="1.007947"/>
  <Type name="1000" class="CT" element="C" mass="12.01078"/>
  <Type name="1001" class="HP" element="H" mass="1.007947"/>
  <Type name="1002" class="CT" element="C" mass="12.01078"/>
  <Type name="1003" class="HC" element="H" mass="1.007947"/>
  <Type name="1004" class="CA" element="C" mass="12.01078"/>
  <Type name="1005" class="CA" element="C" mass="12.01078"/>
  <Type name="1006" class="HA" element="H" mass="1.007947"/>
  <Type name="1007" class="CA" element="C" mass="12.01078"/>
  <Type name="1008" class="HA" element="H" mass="1.007947"/>
  <Type name="1009" class="C" element="C" mass="12.01078"/>
  <Type name="1010" class="OH" element="O" mass="15.99943"/>
  <Type name="1011" class="HO" element="H" mass="1.007947"/>
  <Type name="1012" class="CA" element="C" mass="12.01078"/>
  <Type name="1013" class="HA" element="H" mass="1.007947"/>
  <Type name="1014" class="CA" element="C" mass="12.01078"/>
  <Type name="1015" class="HA" element="H" mass="1.007947"/>
  <Type name="1016" class="C" element="C" mass="12.01078"/>
  <Type name="1017" class="O" element="O" mass="15.99943"/>
  <Type name="1018" class="N3" element="N" mass="14.00672"/>
  <Type name="1019" class="H" element="H" mass="1.007947"/>
  <Type name="1020" class="CT" element="C" mass="12.01078"/>
  <Type name="1021" class="HP" element="H" mass="1.007947"/>
  <Type name="1022" class="CT" element="C" mass="12.01078"/>
  <Type name="1023" class="HC" element="H" mass="1.007947"/>
  <Type name="1024" class="CT" element="C" mass="12.01078"/>
  <Type name="1025" class="HC" element="H" mass="1.007947"/>
  <Type name="1026" class="CT" element="C" mass="12.01078"/>
  <Type name="1027" class="HC" element="H" mass="1.007947"/>
  <Type name="1028" class="C" element="C" mass="12.01078"/>
  <Type name="1029" class="O" element="O" mass="15.99943"/>
  <Type name="1030" class="P" element="P" mass="30.9737622"/>
  <Type name="1031" class="O2" element="O" mass="15.99943"/>
  <Type name="1032" class="O2" element="O" mass="15.99943"/>
  <Type name="1033" class="OS" element="O" mass="15.99943"/>
  <Type name="1034" class="CT" element="C" mass="12.01078"/>
  <Type name="1035" class="H1" element="H" mass="1.007947"/>
  <Type name="1036" class="CT" element="C" mass="12.01078"/>
  <Type name="1037" class="H1" element="H" mass="1.007947"/>
  <Type name="1038" class="OS" element="O" mass="15.99943"/>
  <Type name="1039" class="CT" element="C" mass="12.01078"/>
  <Type name="1040" class="H2" element="H" mass="1.007947"/>
  <Type name="1041" class="N*" element="N" mass="14.00672"/>
  <Type name="1042" class="CK" element="C" mass="12.01078"/>
  <Type name="1043" class="H5" element="H" mass="1.007947"/>
  <Type name="1044" class="NB" element="N" mass="14.00672"/>
  <Type name="1045" class="CB" element="C" mass="12.01078"/>
  <Type name="1046" class="CA" element="C" mass="12.01078"/>
  <Type name="1047" class="N2" element="N" mass="14.00672"/>
  <Type name="1048" class="H" element="H" mass="1.007947"/>
  <Type name="1049" class="NC" element="N" mass="14.00672"/>
  <Type name="1050" class="CQ" element="C" mass="12.01078"/>
  <Type name="1051" class="H5" element="H" mass="1.007947"/>
  <Type name="1052" class="NC" element="N" mass="14.00672"/>
  <Type name="1053" class="CB" element="C" mass="12.01078"/>
  <Type name="1054" class="CT" element="C" mass="12.01078"/>
  <Type name="1055" class="H1" element="H" mass="1.007947"/>
  <Type name="1056" class="CT" element="C" mass="12.01078"/>
  <Type name="1057" class="HC" element="H" mass="1.007947"/>
  <Type name="1058" class="OS" element="O" mass="15.99943"/>
  <Type name="1059" class="P" element="P" mass="30.9737622"/>
  <Type name="1060" class="O2" element="O" mass="15.99943"/>
  <Type name="1061" class="O2" element="O" mass="15.99943"/>
  <Type name="1062" class="OS" element="O" mass="15.99943"/>
  <Type name="1063" class="CT" element="C" mass="12.01078"/>
  <Type name="1064" class="H1" element="H" mass="1.007947"/>
  <Type name="1065" class="CT" element="C" mass="12.01078"/>
  <Type name="1066" class="H1" element="H" mass="1.007947"/>
  <Type name="1067" class="OS" element="O" mass="15.99943"/>
  <Type name="1068" class="CT" element="C" mass="12.01078"/>
  <Type name="1069" class="H2" element="H" mass="1.007947"/>
  <Type name="1070" class="N*" element="N" mass="14.00672"/>
  <Type name="1071" class="CK" element="C" mass="12.01078"/>
  <Type name="1072" class="H5" element="H" mass="1.007947"/>
  <Type name="1073" class="NB" element="N" mass="14.00672"/>
  <Type name="1074" class="CB" element="C" mass="12.01078"/>
  <Type name="1075" class="CA" element="C" mass="12.01078"/>
  <Type name="1076" class="N2" element="N" mass="14.00672"/>
  <Type name="1077" class="H" element="H" mass="1.007947"/>
  <Type name="1078" class="NC" element="N" mass="14.00672"/>
  <Type name="1079" class="CQ" element="C" mass="12.01078"/>
  <Type name="1080" class="H5" element="H" mass="1.007947"/>
  <Type name="1081" class="NC" element="N" mass="14.00672"/>
  <Type name="1082" class="CB" element="C" mass="12.01078"/>
  <Type name="1083" class="CT" element="C" mass="12.01078"/>
  <Type name="1084" class="H1" element="H" mass="1.007947"/>
  <Type name="1085" class="CT" element="C" mass="12.01078"/>
  <Type name="1086" class="HC" element="H" mass="1.007947"/>
  <Type name="1087" class="OH" element="O" mass="15.99943"/>
  <Type name="1088" class="HO" element="H" mass="1.007947"/>
  <Type name="1089" class="HO" element="H" mass="1.007947"/>
  <Type name="1090" class="OH" element="O" mass="15.99943"/>
  <Type name="1091" class="CT" element="C" mass="12.01078"/>
  <Type name="1092" class="H1" element="H" mass="1.007947"/>
  <Type name="1093" class="CT" element="C" mass="12.01078"/>
  <Type name="1094" class="H1" element="H" mass="1.007947"/>
  <Type name="1095" class="OS" element="O" mass="15.99943"/>
  <Type name="1096" class="CT" element="C" mass="12.01078"/>
  <Type name="1097" class="H2" element="H" mass="1.007947"/>
  <Type name="1098" class="N*" element="N" mass="14.00672"/>
  <Type name="1099" class="CK" element="C" mass="12.01078"/>
  <Type name="1100" class="H5" element="H" mass="1.007947"/>
  <Type name="1101" class="NB" element="N" mass="14.00672"/>
  <Type name="1102" class="CB" element="C" mass="12.01078"/>
  <Type name="1103" class="CA" element="C" mass="12.01078"/>
  <Type name="1104" class="N2" element="N" mass="14.00672"/>
  <Type name="1105" class="H" element="H" mass="1.007947"/>
  <Type name="1106" class="NC" element="N" mass="14.00672"/>
  <Type name="1107" class="CQ" element="C" mass="12.01078"/>
  <Type name="1108" class="H5" element="H" mass="1.007947"/>
  <Type name="1109" class="NC" element="N" mass="14.00672"/>
  <Type name="1110" class="CB" element="C" mass="12.01078"/>
  <Type name="1111" class="CT" element="C" mass="12.01078"/>
  <Type name="1112" class="H1" element="H" mass="1.007947"/>
  <Type name="1113" class="CT" element="C" mass="12.01078"/>
  <Type name="1114" class="HC" element="H" mass="1.007947"/>
  <Type name="1115" class="OS" element="O" mass="15.99943"/>
  <Type name="1116" class="HO" element="H" mass="1.007947"/>
  <Type name="1117" class="OH" element="O" mass="15.99943"/>
  <Type name="1118" class="CT" element="C" mass="12.01078"/>
  <Type name="1119" class="H1" element="H" mass="1.007947"/>
  <Type name="1120" class="CT" element="C" mass="12.01078"/>
  <Type name="1121" class="H1" element="H" mass="1.007947"/>
  <Type name="1122" class="OS" element="O" mass="15.99943"/>
  <Type name="1123" class="CT" element="C" mass="12.01078"/>
  <Type name="1124" class="H2" element="H" mass="1.007947"/>
  <Type name="1125" class="N*" element="N" mass="14.00672"/>
  <Type name="1126" class="CK" element="C" mass="12.01078"/>
  <Type name="1127" class="H5" element="H" mass="1.007947"/>
  <Type name="1128" class="NB" element="N" mass="14.00672"/>
  <Type name="1129" class="CB" element="C" mass="12.01078"/>
  <Type name="1130" class="CA" element="C" mass="12.01078"/>
  <Type name="1131" class="N2" element="N" mass="14.00672"/>
  <Type name="1132" class="H" element="H" mass="1.007947"/>
  <Type name="1133" class="NC" element="N" mass="14.00672"/>
  <Type name="1134" class="CQ" element="C" mass="12.01078"/>
  <Type name="1135" class="H5" element="H" mass="1.007947"/>
  <Type name="1136" class="NC" element="N" mass="14.00672"/>
  <Type name="1137" class="CB" element="C" mass="12.01078"/>
  <Type name="1138" class="CT" element="C" mass="12.01078"/>
  <Type name="1139" class="H1" element="H" mass="1.007947"/>
  <Type name="1140" class="CT" element="C" mass="12.01078"/>
  <Type name="1141" class="HC" element="H" mass="1.007947"/>
  <Type name="1142" class="OH" element="O" mass="15.99943"/>
  <Type name="1143" class="HO" element="H" mass="1.007947"/>
  <Type name="1144" class="P" element="P" mass="30.9737622"/>
  <Type name="1145" class="O2" element="O" mass="15.99943"/>
  <Type name="1146" class="O2" element="O" mass="15.99943"/>
  <Type name="1147" class="OS" element="O" mass="15.99943"/>
  <Type name="1148" class="CT" element="C" mass="12.01078"/>
  <Type name="1149" class="H1" element="H" mass="1.007947"/>
  <Type name="1150" class="CT" element="C" mass="12.01078"/>
  <Type name="1151" class="H1" element="H" mass="1.007947"/>
  <Type name="1152" class="OS" element="O" mass="15.99943"/>
  <Type name="1153" class="CT" element="C" mass="12.01078"/>
  <Type name="1154" class="H2" element="H" mass="1.007947"/>
  <Type name="1155" class="N*" element="N" mass="14.00672"/>
  <Type name="1156" class="CM" element="C" mass="12.01078"/>
  <Type name="1157" class="H4" element="H" mass="1.007947"/>
  <Type name="1158" class="CM" element="C" mass="12.01078"/>
  <Type name="1159" class="HA" element="H" mass="1.007947"/>
  <Type name="1160" class="CA" element="C" mass="12.01078"/>
  <Type name="1161" class="N2" element="N" mass="14.00672"/>
  <Type name="1162" class="H" element="H" mass="1.007947"/>
  <Type name="1163" class="NC" element="N" mass="14.00672"/>
  <Type name="1164" class="C" element="C" mass="12.01078"/>
  <Type name="1165" class="O" element="O" mass="15.99943"/>
  <Type name="1166" class="CT" element="C" mass="12.01078"/>
  <Type name="1167" class="H1" element="H" mass="1.007947"/>
  <Type name="1168" class="CT" element="C" mass="12.01078"/>
  <Type name="1169" class="HC" element="H" mass="1.007947"/>
  <Type name="1170" class="OS" element="O" mass="15.99943"/>
  <Type name="1171" class="P" element="P" mass="30.9737622"/>
  <Type name="1172" class="O2" element="O" mass="15.99943"/>
  <Type name="1173" class="O2" element="O" mass="15.99943"/>
  <Type name="1174" class="OS" element="O" mass="15.99943"/>
  <Type name="1175" class="CT" element="C" mass="12.01078"/>
  <Type name="1176" class="H1" element="H" mass="1.007947"/>
  <Type name="1177" class="CT" element="C" mass="12.01078"/>
  <Type name="1178" class="H1" element="H" mass="1.007947"/>
  <Type name="1179" class="OS" element="O" mass="15.99943"/>
  <Type name="1180" class="CT" element="C" mass="12.01078"/>
  <Type name="1181" class="H2" element="H" mass="1.007947"/>
  <Type name="1182" class="N*" element="N" mass="14.00672"/>
  <Type name="1183" class="CM" element="C" mass="12.01078"/>
  <Type name="1184" class="H4" element="H" mass="1.007947"/>
  <Type name="1185" class="CM" element="C" mass="12.01078"/>
  <Type name="1186" class="HA" element="H" mass="1.007947"/>
  <Type name="1187" class="CA" element="C" mass="12.01078"/>
  <Type name="1188" class="N2" element="N" mass="14.00672"/>
  <Type name="1189" class="H" element="H" mass="1.007947"/>
  <Type name="1190" class="NC" element="N" mass="14.00672"/>
  <Type name="1191" class="C" element="C" mass="12.01078"/>
  <Type name="1192" class="O" element="O" mass="15.99943"/>
  <Type name="1193" class="CT" element="C" mass="12.01078"/>
  <Type name="1194" class="H1" element="H" mass="1.007947"/>
  <Type name="1195" class="CT" element="C" mass="12.01078"/>
  <Type name="1196" class="HC" element="H" mass="1.007947"/>
  <Type name="1197" class="OH" element="O" mass="15.99943"/>
  <Type name="1198" class="HO" element="H" mass="1.007947"/>
  <Type name="1199" class="HO" element="H" mass="1.007947"/>
  <Type name="1200" class="OH" element="O" mass="15.99943"/>
  <Type name="1201" class="CT" element="C" mass="12.01078"/>
  <Type name="1202" class="H1" element="H" mass="1.007947"/>
  <Type name="1203" class="CT" element="C" mass="12.01078"/>
  <Type name="1204" class="H1" element="H" mass="1.007947"/>
  <Type name="1205" class="OS" element="O" mass="15.99943"/>
  <Type name="1206" class="CT" element="C" mass="12.01078"/>
  <Type name="1207" class="H2" element="H" mass="1.007947"/>
  <Type name="1208" class="N*" element="N" mass="14.00672"/>
  <Type name="1209" class="CM" element="C" mass="12.01078"/>
  <Type name="1210" class="H4" element="H" mass="1.007947"/>
  <Type name="1211" class="CM" element="C" mass="12.01078"/>
  <Type name="1212" class="HA" element="H" mass="1.007947"/>
  <Type name="1213" class="CA" element="C" mass="12.01078"/>
  <Type name="1214" class="N2" element="N" mass="14.00672"/>
  <Type name="1215" class="H" element="H" mass="1.007947"/>
  <Type name="1216" class="NC" element="N" mass="14.00672"/>
  <Type name="1217" class="C" element="C" mass="12.01078"/>
  <Type name="1218" class="O" element="O" mass="15.99943"/>
  <Type name="1219" class="CT" element="C" mass="12.01078"/>
  <Type name="1220" class="H1" element="H" mass="1.007947"/>
  <Type name="1221" class="CT" element="C" mass="12.01078"/>
  <Type name="1222" class="HC" element="H" mass="1.007947"/>
  <Type name="1223" class="OS" element="O" mass="15.99943"/>
  <Type name="1224" class="HO" element="H" mass="1.007947"/>
  <Type name="1225" class="OH" element="O" mass="15.99943"/>
  <Type name="1226" class="CT" element="C" mass="12.01078"/>
  <Type name="1227" class="H1" element="H" mass="1.007947"/>
  <Type name="1228" class="CT" element="C" mass="12.01078"/>
  <Type name="1229" class="H1" element="H" mass="1.007947"/>
  <Type name="1230" class="OS" element="O" mass="15.99943"/>
  <Type name="1231" class="CT" element="C" mass="12.01078"/>
  <Type name="1232" class="H2" element="H" mass="1.007947"/>
  <Type name="1233" class="N*" element="N" mass="14.00672"/>
  <Type name="1234" class="CM" element="C" mass="12.01078"/>
  <Type name="1235" class="H4" element="H" mass="1.007947"/>
  <Type name="1236" class="CM" element="C" mass="12.01078"/>
  <Type name="1237" class="HA" element="H" mass="1.007947"/>
  <Type name="1238" class="CA" element="C" mass="12.01078"/>
  <Type name="1239" class="N2" element="N" mass="14.00672"/>
  <Type name="1240" class="H" element="H" mass="1.007947"/>
  <Type name="1241" class="NC" element="N" mass="14.00672"/>
  <Type name="1242" class="C" element="C" mass="12.01078"/>
  <Type name="1243" class="O" element="O" mass="15.99943"/>
  <Type name="1244" class="CT" element="C" mass="12.01078"/>
  <Type name="1245" class="H1" element="H" mass="1.007947"/>
  <Type name="1246" class="CT" element="C" mass="12.01078"/>
  <Type name="1247" class="HC" element="H" mass="1.007947"/>
  <Type name="1248" class="OH" element="O" mass="15.99943"/>
  <Type name="1249" class="HO" element="H" mass="1.007947"/>
  <Type name="1250" class="P" element="P" mass="30.9737622"/>
  <Type name="1251" class="O2" element="O" mass="15.99943"/>
  <Type name="1252" class="O2" element="O" mass="15.99943"/>
  <Type name="1253" class="OS" element="O" mass="15.99943"/>
  <Type name="1254" class="CT" element="C" mass="12.01078"/>
  <Type name="1255" class="H1" element="H" mass="1.007947"/>
  <Type name="1256" class="CT" element="C" mass="12.01078"/>
  <Type name="1257" class="H1" element="H" mass="1.007947"/>
  <Type name="1258" class="OS" element="O" mass="15.99943"/>
  <Type name="1259" class="CT" element="C" mass="12.01078"/>
  <Type name="1260" class="H2" element="H" mass="1.007947"/>
  <Type name="1261" class="N*" element="N" mass="14.00672"/>
  <Type name="1262" class="CK" element="C" mass="12.01078"/>
  <Type name="1263" class="H5" element="H" mass="1.007947"/>
  <Type name="1264" class="NB" element="N" mass="14.00672"/>
  <Type name="1265" class="CB" element="C" mass="12.01078"/>
  <Type name="1266" class="C" element="C" mass="12.01078"/>
  <Type name="1267" class="O" element="O" mass="15.99943"/>
  <Type name="1268" class="NA" element="N" mass="14.00672"/>
  <Type name="1269" class="H" element="H" mass="1.007947"/>
  <Type name="1270" class="CA" element="C" mass="12.01078"/>
  <Type name="1271" class="N2" element="N" mass="14.00672"/>
  <Type name="1272" class="H" element="H" mass="1.007947"/>
  <Type name="1273" class="NC" element="N" mass="14.00672"/>
  <Type name="1274" class="CB" element="C" mass="12.01078"/>
  <Type name="1275" class="CT" element="C" mass="12.01078"/>
  <Type name="1276" class="H1" element="H" mass="1.007947"/>
  <Type name="1277" class="CT" element="C" mass="12.01078"/>
  <Type name="1278" class="HC" element="H" mass="1.007947"/>
  <Type name="1279" class="OS" element="O" mass="15.99943"/>
  <Type name="1280" class="P" element="P" mass="30.9737622"/>
  <Type name="1281" class="O2" element="O" mass="15.99943"/>
  <Type name="1282" class="O2" element="O" mass="15.99943"/>
  <Type name="1283" class="OS" element="O" mass="15.99943"/>
  <Type name="1284" class="CT" element="C" mass="12.01078"/>
  <Type name="1285" class="H1" element="H" mass="1.007947"/>
  <Type name="1286" class="CT" element="C" mass="12.01078"/>
  <Type name="1287" class="H1" element="H" mass="1.007947"/>
  <Type name="1288" class="OS" element="O" mass="15.99943"/>
  <Type name="1289" class="CT" element="C" mass="12.01078"/>
  <Type name="1290" class="H2" element="H" mass="1.007947"/>
  <Type name="1291" class="N*" element="N" mass="14.00672"/>
  <Type name="1292" class="CK" element="C" mass="12.01078"/>
  <Type name="1293" class="H5" element="H" mass="1.007947"/>
  <Type name="1294" class="NB" element="N" mass="14.00672"/>
  <Type name="1295" class="CB" element="C" mass="12.01078"/>
  <Type name="1296" class="C" element="C" mass="12.01078"/>
  <Type name="1297" class="O" element="O" mass="15.99943"/>
  <Type name="1298" class="NA" element="N" mass="14.00672"/>
  <Type name="1299" class="H" element="H" mass="1.007947"/>
  <Type name="1300" class="CA" element="C" mass="12.01078"/>
  <Type name="1301" class="N2" element="N" mass="14.00672"/>
  <Type name="1302" class="H" element="H" mass="1.007947"/>
  <Type name="1303" class="NC" element="N" mass="14.00672"/>
  <Type name="1304" class="CB" element="C" mass="12.01078"/>
  <Type name="1305" class="CT" element="C" mass="12.01078"/>
  <Type name="1306" class="H1" element="H" mass="1.007947"/>
  <Type name="1307" class="CT" element="C" mass="12.01078"/>
  <Type name="1308" class="HC" element="H" mass="1.007947"/>
  <Type name="1309" class="OH" element="O" mass="15.99943"/>
  <Type name="1310" class="HO" element="H" mass="1.007947"/>
  <Type name="1311" class="HO" element="H" mass="1.007947"/>
  <Type name="1312" class="OH" element="O" mass="15.99943"/>
  <Type name="1313" class="CT" element="C" mass="12.01078"/>
  <Type name="1314" class="H1" element="H" mass="1.007947"/>
  <Type name="1315" class="CT" element="C" mass="12.01078"/>
  <Type name="1316" class="H1" element="H" mass="1.007947"/>
  <Type name="1317" class="OS" element="O" mass="15.99943"/>
  <Type name="1318" class="CT" element="C" mass="12.01078"/>
  <Type name="1319" class="H2" element="H" mass="1.007947"/>
  <Type name="1320" class="N*" element="N" mass="14.00672"/>
  <Type name="1321" class="CK" element="C" mass="12.01078"/>
  <Type name="1322" class="H5" element="H" mass="1.007947"/>
  <Type name="1323" class="NB" element="N" mass="14.00672"/>
  <Type name="1324" class="CB" element="C" mass="12.01078"/>
  <Type name="1325" class="C" element="C" mass="12.01078"/>
  <Type name="1326" class="O" element="O" mass="15.99943"/>
  <Type name="1327" class="NA" element="N" mass="14.00672"/>
  <Type name="1328" class="H" element="H" mass="1.007947"/>
  <Type name="1329" class="CA" element="C" mass="12.01078"/>
  <Type name="1330" class="N2" element="N" mass="14.00672"/>
  <Type name="1331" class="H" element="H" mass="1.007947"/>
  <Type name="1332" class="NC" element="N" mass="14.00672"/>
  <Type name="1333" class="CB" element="C" mass="12.01078"/>
  <Type name="1334" class="CT" element="C" mass="12.01078"/>
  <Type name="1335" class="H1" element="H" mass="1.007947"/>
  <Type name="1336" class="CT" element="C" mass="12.01078"/>
  <Type name="1337" class="HC" element="H" mass="1.007947"/>
  <Type name="1338" class="OS" element="O" mass="15.99943"/>
  <Type name="1339" class="HO" element="H" mass="1.007947"/>
  <Type name="1340" class="OH" element="O" mass="15.99943"/>
  <Type name="1341" class="CT" element="C" mass="12.01078"/>
  <Type name="1342" class="H1" element="H" mass="1.007947"/>
  <Type name="1343" class="CT" element="C" mass="12.01078"/>
  <Type name="1344" class="H1" element="H" mass="1.007947"/>
  <Type name="1345" class="OS" element="O" mass="15.99943"/>
  <Type name="1346" class="CT" element="C" mass="12.01078"/>
  <Type name="1347" class="H2" element="H" mass="1.007947"/>
  <Type name="1348" class="N*" element="N" mass="14.00672"/>
  <Type name="1349" class="CK" element="C" mass="12.01078"/>
  <Type name="1350" class="H5" element="H" mass="1.007947"/>
  <Type name="1351" class="NB" element="N" mass="14.00672"/>
  <Type name="1352" class="CB" element="C" mass="12.01078"/>
  <Type name="1353" class="C" element="C" mass="12.01078"/>
  <Type name="1354" class="O" element="O" mass="15.99943"/>
  <Type name="1355" class="NA" element="N" mass="14.00672"/>
  <Type name="1356" class="H" element="H" mass="1.007947"/>
  <Type name="1357" class="CA" element="C" mass="12.01078"/>
  <Type name="1358" class="N2" element="N" mass="14.00672"/>
  <Type name="1359" class="H" element="H" mass="1.007947"/>
  <Type name="1360" class="NC" element="N" mass="14.00672"/>
  <Type name="1361" class="CB" element="C" mass="12.01078"/>
  <Type name="1362" class="CT" element="C" mass="12.01078"/>
  <Type name="1363" class="H1" element="H" mass="1.007947"/>
  <Type name="1364" class="CT" element="C" mass="12.01078"/>
  <Type name="1365" class="HC" element="H" mass="1.007947"/>
  <Type name="1366" class="OH" element="O" mass="15.99943"/>
  <Type name="1367" class="HO" element="H" mass="1.007947"/>
  <Type name="1368" class="P" element="P" mass="30.9737622"/>
  <Type name="1369" class="O2" element="O" mass="15.99943"/>
  <Type name="1370" class="O2" element="O" mass="15.99943"/>
  <Type name="1371" class="OS" element="O" mass="15.99943"/>
  <Type name="1372" class="CT" element="C" mass="12.01078"/>
  <Type name="1373" class="H1" element="H" mass="1.007947"/>
  <Type name="1374" class="CT" element="C" mass="12.01078"/>
  <Type name="1375" class="H1" element="H" mass="1.007947"/>
  <Type name="1376" class="OS" element="O" mass="15.99943"/>
  <Type name="1377" class="CT" element="C" mass="12.01078"/>
  <Type name="1378" class="H2" element="H" mass="1.007947"/>
  <Type name="1379" class="N*" element="N" mass="14.00672"/>
  <Type name="1380" class="CM" element="C" mass="12.01078"/>
  <Type name="1381" class="H4" element="H" mass="1.007947"/>
  <Type name="1382" class="CM" element="C" mass="12.01078"/>
  <Type name="1383" class="CT" element="C" mass="12.01078"/>
  <Type name="1384" class="HC" element="H" mass="1.007947"/>
  <Type name="1385" class="C" element="C" mass="12.01078"/>
  <Type name="1386" class="O" element="O" mass="15.99943"/>
  <Type name="1387" class="NA" element="N" mass="14.00672"/>
  <Type name="1388" class="H" element="H" mass="1.007947"/>
  <Type name="1389" class="C" element="C" mass="12.01078"/>
  <Type name="1390" class="O" element="O" mass="15.99943"/>
  <Type name="1391" class="CT" element="C" mass="12.01078"/>
  <Type name="1392" class="H1" element="H" mass="1.007947"/>
  <Type name="1393" class="CT" element="C" mass="12.01078"/>
  <Type name="1394" class="HC" element="H" mass="1.007947"/>
  <Type name="1395" class="OS" element="O" mass="15.99943"/>
  <Type name="1396" class="P" element="P" mass="30.9737622"/>
  <Type name="1397" class="O2" element="O" mass="15.99943"/>
  <Type name="1398" class="O2" element="O" mass="15.99943"/>
  <Type name="1399" class="OS" element="O" mass="15.99943"/>
  <Type name="1400" class="CT" element="C" mass="12.01078"/>
  <Type name="1401" class="H1" element="H" mass="1.007947"/>
  <Type name="1402" class="CT" element="C" mass="12.01078"/>
  <Type name="1403" class="H1" element="H" mass="1.007947"/>
  <Type name="1404" class="OS" element="O" mass="15.99943"/>
  <Type name="1405" class="CT" element="C" mass="12.01078"/>
  <Type name="1406" class="H2" element="H" mass="1.007947"/>
  <Type name="1407" class="N*" element="N" mass="14.00672"/>
  <Type name="1408" class="CM" element="C" mass="12.01078"/>
  <Type name="1409" class="H4" element="H" mass="1.007947"/>
  <Type name="1410" class="CM" element="C" mass="12.01078"/>
  <Type name="1411" class="CT" element="C" mass="12.01078"/>
  <Type name="1412" class="HC" element="H" mass="1.007947"/>
  <Type name="1413" class="C" element="C" mass="12.01078"/>
  <Type name="1414" class="O" element="O" mass="15.99943"/>
  <Type name="1415" class="NA" element="N" mass="14.00672"/>
  <Type name="1416" class="H" element="H" mass="1.007947"/>
  <Type name="1417" class="C" element="C" mass="12.01078"/>
  <Type name="1418" class="O" element="O" mass="15.99943"/>
  <Type name="1419" class="CT" element="C" mass="12.01078"/>
  <Type name="1420" class="H1" element="H" mass="1.007947"/>
  <Type name="1421" class="CT" element="C" mass="12.01078"/>
  <Type name="1422" class="HC" element="H" mass="1.007947"/>
  <Type name="1423" class="OH" element="O" mass="15.99943"/>
  <Type name="1424" class="HO" element="H" mass="1.007947"/>
  <Type name="1425" class="HO" element="H" mass="1.007947"/>
  <Type name="1426" class="OH" element="O" mass="15.99943"/>
  <Type name="1427" class="CT" element="C" mass="12.01078"/>
  <Type name="1428" class="H1" element="H" mass="1.007947"/>
  <Type name="1429" class="CT" element="C" mass="12.01078"/>
  <Type name="1430" class="H1" element="H" mass="1.007947"/>
  <Type name="1431" class="OS" element="O" mass="15.99943"/>
  <Type name="1432" class="CT" element="C" mass="12.01078"/>
  <Type name="1433" class="H2" element="H" mass="1.007947"/>
  <Type name="1434" class="N*" element="N" mass="14.00672"/>
  <Type name="1435" class="CM" element="C" mass="12.01078"/>
  <Type name="1436" class="H4" element="H" mass="1.007947"/>
  <Type name="1437" class="CM" element="C" mass="12.01078"/>
  <Type name="1438" class="CT" element="C" mass="12.01078"/>
  <Type name="1439" class="HC" element="H" mass="1.007947"/>
  <Type name="1440" class="C" element="C" mass="12.01078"/>
  <Type name="1441" class="O" element="O" mass="15.99943"/>
  <Type name="1442" class="NA" element="N" mass="14.00672"/>
  <Type name="1443" class="H" element="H" mass="1.007947"/>
  <Type name="1444" class="C" element="C" mass="12.01078"/>
  <Type name="1445" class="O" element="O" mass="15.99943"/>
  <Type name="1446" class="CT" element="C" mass="12.01078"/>
  <Type name="1447" class="H1" element="H" mass="1.007947"/>
  <Type name="1448" class="CT" element="C" mass="12.01078"/>
  <Type name="1449" class="HC" element="H" mass="1.007947"/>
  <Type name="1450" class="OS" element="O" mass="15.99943"/>
  <Type name="1451" class="HO" element="H" mass="1.007947"/>
  <Type name="1452" class="OH" element="O" mass="15.99943"/>
  <Type name="1453" class="CT" element="C" mass="12.01078"/>
  <Type name="1454" class="H1" element="H" mass="1.007947"/>
  <Type name="1455" class="CT" element="C" mass="12.01078"/>
  <Type name="1456" class="H1" element="H" mass="1.007947"/>
  <Type name="1457" class="OS" element="O" mass="15.99943"/>
  <Type name="1458" class="CT" element="C" mass="12.01078"/>
  <Type name="1459" class="H2" element="H" mass="1.007947"/>
  <Type name="1460" class="N*" element="N" mass="14.00672"/>
  <Type name="1461" class="CM" element="C" mass="12.01078"/>
  <Type name="1462" class="H4" element="H" mass="1.007947"/>
  <Type name="1463" class="CM" element="C" mass="12.01078"/>
  <Type name="1464" class="CT" element="C" mass="12.01078"/>
  <Type name="1465" class="HC" element="H" mass="1.007947"/>
  <Type name="1466" class="C" element="C" mass="12.01078"/>
  <Type name="1467" class="O" element="O" mass="15.99943"/>
  <Type name="1468" class="NA" element="N" mass="14.00672"/>
  <Type name="1469" class="H" element="H" mass="1.007947"/>
  <Type name="1470" class="C" element="C" mass="12.01078"/>
  <Type name="1471" class="O" element="O" mass="15.99943"/>
  <Type name="1472" class="CT" element="C" mass="12.01078"/>
  <Type name="1473" class="H1" element="H" mass="1.007947"/>
  <Type name="1474" class="CT" element="C" mass="12.01078"/>
  <Type name="1475" class="HC" element="H" mass="1.007947"/>
  <Type name="1476" class="OH" element="O" mass="15.99943"/>
  <Type name="1477" class="HO" element="H" mass="1.007947"/>
  <Type name="1478" class="P" element="P" mass="30.9737622"/>
  <Type name="1479" class="O2" element="O" mass="15.99943"/>
  <Type name="1480" class="O2" element="O" mass="15.99943"/>
  <Type name="1481" class="OS" element="O" mass="15.99943"/>
  <Type name="1482" class="CT" element="C" mass="12.01078"/>
  <Type name="1483" class="H1" element="H" mass="1.007947"/>
  <Type name="1484" class="CT" element="C" mass="12.01078"/>
  <Type name="1485" class="H1" element="H" mass="1.007947"/>
  <Type name="1486" class="OS" element="O" mass="15.99943"/>
  <Type name="1487" class="CT" element="C" mass="12.01078"/>
  <Type name="1488" class="H2" element="H" mass="1.007947"/>
  <Type name="1489" class="N*" element="N" mass="14.00672"/>
  <Type name="1490" class="CK" element="C" mass="12.01078"/>
  <Type name="1491" class="H5" element="H" mass="1.007947"/>
  <Type name="1492" class="NB" element="N" mass="14.00672"/>
  <Type name="1493" class="CB" element="C" mass="12.01078"/>
  <Type name="1494" class="CA" element="C" mass="12.01078"/>
  <Type name="1495" class="N2" element="N" mass="14.00672"/>
  <Type name="1496" class="H" element="H" mass="1.007947"/>
  <Type name="1497" class="NC" element="N" mass="14.00672"/>
  <Type name="1498" class="CQ" element="C" mass="12.01078"/>
  <Type name="1499" class="H5" element="H" mass="1.007947"/>
  <Type name="1500" class="NC" element="N" mass="14.00672"/>
  <Type name="1501" class="CB" element="C" mass="12.01078"/>
  <Type name="1502" class="CT" element="C" mass="12.01078"/>
  <Type name="1503" class="H1" element="H" mass="1.007947"/>
  <Type name="1504" class="CT" element="C" mass="12.01078"/>
  <Type name="1505" class="H1" element="H" mass="1.007947"/>
  <Type name="1506" class="OH" element="O" mass="15.99943"/>
  <Type name="1507" class="HO" element="H" mass="1.007947"/>
  <Type name="1508" class="OS" element="O" mass="15.99943"/>
  <Type name="1509" class="P" element="P" mass="30.9737622"/>
  <Type name="1510" class="O2" element="O" mass="15.99943"/>
  <Type name="1511" class="O2" element="O" mass="15.99943"/>
  <Type name="1512" class="OS" element="O" mass="15.99943"/>
  <Type name="1513" class="CT" element="C" mass="12.01078"/>
  <Type name="1514" class="H1" element="H" mass="1.007947"/>
  <Type name="1515" class="CT" element="C" mass="12.01078"/>
  <Type name="1516" class="H1" element="H" mass="1.007947"/>
  <Type name="1517" class="OS" element="O" mass="15.99943"/>
  <Type name="1518" class="CT" element="C" mass="12.01078"/>
  <Type name="1519" class="H2" element="H" mass="1.007947"/>
  <Type name="1520" class="N*" element="N" mass="14.00672"/>
  <Type name="1521" class="CK" element="C" mass="12.01078"/>
  <Type name="1522" class="H5" element="H" mass="1.007947"/>
  <Type name="1523" class="NB" element="N" mass="14.00672"/>
  <Type name="1524" class="CB" element="C" mass="12.01078"/>
  <Type name="1525" class="CA" element="C" mass="12.01078"/>
  <Type name="1526" class="N2" element="N" mass="14.00672"/>
  <Type name="1527" class="H" element="H" mass="1.007947"/>
  <Type name="1528" class="NC" element="N" mass="14.00672"/>
  <Type name="1529" class="CQ" element="C" mass="12.01078"/>
  <Type name="1530" class="H5" element="H" mass="1.007947"/>
  <Type name="1531" class="NC" element="N" mass="14.00672"/>
  <Type name="1532" class="CB" element="C" mass="12.01078"/>
  <Type name="1533" class="CT" element="C" mass="12.01078"/>
  <Type name="1534" class="H1" element="H" mass="1.007947"/>
  <Type name="1535" class="CT" element="C" mass="12.01078"/>
  <Type name="1536" class="H1" element="H" mass="1.007947"/>
  <Type name="1537" class="OH" element="O" mass="15.99943"/>
  <Type name="1538" class="HO" element="H" mass="1.007947"/>
  <Type name="1539" class="OH" element="O" mass="15.99943"/>
  <Type name="1540" class="HO" element="H" mass="1.007947"/>
  <Type name="1541" class="HO" element="H" mass="1.007947"/>
  <Type name="1542" class="OH" element="O" mass="15.99943"/>
  <Type name="1543" class="CT" element="C" mass="12.01078"/>
  <Type name="1544" class="H1" element="H" mass="1.007947"/>
  <Type name="1545" class="CT" element="C" mass="12.01078"/>
  <Type name="1546" class="H1" element="H" mass="1.007947"/>
  <Type name="1547" class="OS" element="O" mass="15.99943"/>
  <Type name="1548" class="CT" element="C" mass="12.01078"/>
  <Type name="1549" class="H2" element="H" mass="1.007947"/>
  <Type name="1550" class="N*" element="N" mass="14.00672"/>
  <Type name="1551" class="CK" element="C" mass="12.01078"/>
  <Type name="1552" class="H5" element="H" mass="1.007947"/>
  <Type name="1553" class="NB" element="N" mass="14.00672"/>
  <Type name="1554" class="CB" element="C" mass="12.01078"/>
  <Type name="1555" class="CA" element="C" mass="12.01078"/>
  <Type name="1556" class="N2" element="N" mass="14.00672"/>
  <Type name="1557" class="H" element="H" mass="1.007947"/>
  <Type name="1558" class="NC" element="N" mass="14.00672"/>
  <Type name="1559" class="CQ" element="C" mass="12.01078"/>
  <Type name="1560" class="H5" element="H" mass="1.007947"/>
  <Type name="1561" class="NC" element="N" mass="14.00672"/>
  <Type name="1562" class="CB" element="C" mass="12.01078"/>
  <Type name="1563" class="CT" element="C" mass="12.01078"/>
  <Type name="1564" class="H1" element="H" mass="1.007947"/>
  <Type name="1565" class="CT" element="C" mass="12.01078"/>
  <Type name="1566" class="H1" element="H" mass="1.007947"/>
  <Type name="1567" class="OH" element="O" mass="15.99943"/>
  <Type name="1568" class="HO" element="H" mass="1.007947"/>
  <Type name="1569" class="OS" element="O" mass="15.99943"/>
  <Type name="1570" class="HO" element="H" mass="1.007947"/>
  <Type name="1571" class="OH" element="O" mass="15.99943"/>
  <Type name="1572" class="CT" element="C" mass="12.01078"/>
  <Type name="1573" class="H1" element="H" mass="1.007947"/>
  <Type name="1574" class="CT" element="C" mass="12.01078"/>
  <Type name="1575" class="H1" element="H" mass="1.007947"/>
  <Type name="1576" class="OS" element="O" mass="15.99943"/>
  <Type name="1577" class="CT" element="C" mass="12.01078"/>
  <Type name="1578" class="H2" element="H" mass="1.007947"/>
  <Type name="1579" class="N*" element="N" mass="14.00672"/>
  <Type name="1580" class="CK" element="C" mass="12.01078"/>
  <Type name="1581" class="H5" element="H" mass="1.007947"/>
  <Type name="1582" class="NB" element="N" mass="14.00672"/>
  <Type name="1583" class="CB" element="C" mass="12.01078"/>
  <Type name="1584" class="CA" element="C" mass="12.01078"/>
  <Type name="1585" class="N2" element="N" mass="14.00672"/>
  <Type name="1586" class="H" element="H" mass="1.007947"/>
  <Type name="1587" class="NC" element="N" mass="14.00672"/>
  <Type name="1588" class="CQ" element="C" mass="12.01078"/>
  <Type name="1589" class="H5" element="H" mass="1.007947"/>
  <Type name="1590" class="NC" element="N" mass="14.00672"/>
  <Type name="1591" class="CB" element="C" mass="12.01078"/>
  <Type name="1592" class="CT" element="C" mass="12.01078"/>
  <Type name="1593" class="H1" element="H" mass="1.007947"/>
  <Type name="1594" class="CT" element="C" mass="12.01078"/>
  <Type name="1595" class="H1" element="H" mass="1.007947"/>
  <Type name="1596" class="OH" element="O" mass="15.99943"/>
  <Type name="1597" class="HO" element="H" mass="1.007947"/>
  <Type name="1598" class="OH" element="O" mass="15.99943"/>
  <Type name="1599" class="HO" element="H" mass="1.007947"/>
  <Type name="1600" class="P" element="P" mass="30.9737622"/>
  <Type name="1601" class="O2" element="O" mass="15.99943"/>
  <Type name="1602" class="O2" element="O" mass="15.99943"/>
  <Type name="1603" class="OS" element="O" mass="15.99943"/>
  <Type name="1604" class="CT" element="C" mass="12.01078"/>
  <Type name="1605" class="H1" element="H" mass="1.007947"/>
  <Type name="1606" class="CT" element="C" mass="12.01078"/>
  <Type name="1607" class="H1" element="H" mass="1.007947"/>
  <Type name="1608" class="OS" element="O" mass="15.99943"/>
  <Type name="1609" class="CT" element="C" mass="12.01078"/>
  <Type name="1610" class="H2" element="H" mass="1.007947"/>
  <Type name="1611" class="N*" element="N" mass="14.00672"/>
  <Type name="1612" class="CM" element="C" mass="12.01078"/>
  <Type name="1613" class="H4" element="H" mass="1.007947"/>
  <Type name="1614" class="CM" element="C" mass="12.01078"/>
  <Type name="1615" class="HA" element="H" mass="1.007947"/>
  <Type name="1616" class="CA" element="C" mass="12.01078"/>
  <Type name="1617" class="N2" element="N" mass="14.00672"/>
  <Type name="1618" class="H" element="H" mass="1.007947"/>
  <Type name="1619" class="NC" element="N" mass="14.00672"/>
  <Type name="1620" class="C" element="C" mass="12.01078"/>
  <Type name="1621" class="O" element="O" mass="15.99943"/>
  <Type name="1622" class="CT" element="C" mass="12.01078"/>
  <Type name="1623" class="H1" element="H" mass="1.007947"/>
  <Type name="1624" class="CT" element="C" mass="12.01078"/>
  <Type name="1625" class="H1" element="H" mass="1.007947"/>
  <Type name="1626" class="OH" element="O" mass="15.99943"/>
  <Type name="1627" class="HO" element="H" mass="1.007947"/>
  <Type name="1628" class="OS" element="O" mass="15.99943"/>
  <Type name="1629" class="P" element="P" mass="30.9737622"/>
  <Type name="1630" class="O2" element="O" mass="15.99943"/>
  <Type name="1631" class="O2" element="O" mass="15.99943"/>
  <Type name="1632" class="OS" element="O" mass="15.99943"/>
  <Type name="1633" class="CT" element="C" mass="12.01078"/>
  <Type name="1634" class="H1" element="H" mass="1.007947"/>
  <Type name="1635" class="CT" element="C" mass="12.01078"/>
  <Type name="1636" class="H1" element="H" mass="1.007947"/>
  <Type name="1637" class="OS" element="O" mass="15.99943"/>
  <Type name="1638" class="CT" element="C" mass="12.01078"/>
  <Type name="1639" class="H2" element="H" mass="1.007947"/>
  <Type name="1640" class="N*" element="N" mass="14.00672"/>
  <Type name="1641" class="CM" element="C" mass="12.01078"/>
  <Type name="1642" class="H4" element="H" mass="1.007947"/>
  <Type name="1643" class="CM" element="C" mass="12.01078"/>
  <Type name="1644" class="HA" element="H" mass="1.007947"/>
  <Type name="1645" class="CA" element="C" mass="12.01078"/>
  <Type name="1646" class="N2" element="N" mass="14.00672"/>
  <Type name="1647" class="H" element="H" mass="1.007947"/>
  <Type name="1648" class="NC" element="N" mass="14.00672"/>
  <Type name="1649" class="C" element="C" mass="12.01078"/>
  <Type name="1650" class="O" element="O" mass="15.99943"/>
  <Type name="1651" class="CT" element="C" mass="12.01078"/>
  <Type name="1652" class="H1" element="H" mass="1.007947"/>
  <Type name="1653" class="CT" element="C" mass="12.01078"/>
  <Type name="1654" class="H1" element="H" mass="1.007947"/>
  <Type name="1655" class="OH" element="O" mass="15.99943"/>
  <Type name="1656" class="HO" element="H" mass="1.007947"/>
  <Type name="1657" class="OH" element="O" mass="15.99943"/>
  <Type name="1658" class="HO" element="H" mass="1.007947"/>
  <Type name="1659" class="HO" element="H" mass="1.007947"/>
  <Type name="1660" class="OH" element="O" mass="15.99943"/>
  <Type name="1661" class="CT" element="C" mass="12.01078"/>
  <Type name="1662" class="H1" element="H" mass="1.007947"/>
  <Type name="1663" class="CT" element="C" mass="12.01078"/>
  <Type name="1664" class="H1" element="H" mass="1.007947"/>
  <Type name="1665" class="OS" element="O" mass="15.99943"/>
  <Type name="1666" class="CT" element="C" mass="12.01078"/>
  <Type name="1667" class="H2" element="H" mass="1.007947"/>
  <Type name="1668" class="N*" element="N" mass="14.00672"/>
  <Type name="1669" class="CM" element="C" mass="12.01078"/>
  <Type name="1670" class="H4" element="H" mass="1.007947"/>
  <Type name="1671" class="CM" element="C" mass="12.01078"/>
  <Type name="1672" class="HA" element="H" mass="1.007947"/>
  <Type name="1673" class="CA" element="C" mass="12.01078"/>
  <Type name="1674" class="N2" element="N" mass="14.00672"/>
  <Type name="1675" class="H" element="H" mass="1.007947"/>
  <Type name="1676" class="NC" element="N" mass="14.00672"/>
  <Type name="1677" class="C" element="C" mass="12.01078"/>
  <Type name="1678" class="O" element="O" mass="15.99943"/>
  <Type name="1679" class="CT" element="C" mass="12.01078"/>
  <Type name="1680" class="H1" element="H" mass="1.007947"/>
  <Type name="1681" class="CT" element="C" mass="12.01078"/>
  <Type name="1682" class="H1" element="H" mass="1.007947"/>
  <Type name="1683" class="OH" element="O" mass="15.99943"/>
  <Type name="1684" class="HO" element="H" mass="1.007947"/>
  <Type name="1685" class="OS" element="O" mass="15.99943"/>
  <Type name="1686" class="HO" element="H" mass="1.007947"/>
  <Type name="1687" class="OH" element="O" mass="15.99943"/>
  <Type name="1688" class="CT" element="C" mass="12.01078"/>
  <Type name="1689" class="H1" element="H" mass="1.007947"/>
  <Type name="1690" class="CT" element="C" mass="12.01078"/>
  <Type name="1691" class="H1" element="H" mass="1.007947"/>
  <Type name="1692" class="OS" element="O" mass="15.99943"/>
  <Type name="1693" class="CT" element="C" mass="12.01078"/>
  <Type name="1694" class="H2" element="H" mass="1.007947"/>
  <Type name="1695" class="N*" element="N" mass="14.00672"/>
  <Type name="1696" class="CM" element="C" mass="12.01078"/>
  <Type name="1697" class="H4" element="H" mass="1.007947"/>
  <Type name="1698" class="CM" element="C" mass="12.01078"/>
  <Type name="1699" class="HA" element="H" mass="1.007947"/>
  <Type name="1700" class="CA" element="C" mass="12.01078"/>
  <Type name="1701" class="N2" element="N" mass="14.00672"/>
  <Type name="1702" class="H" element="H" mass="1.007947"/>
  <Type name="1703" class="NC" element="N" mass="14.00672"/>
  <Type name="1704" class="C" element="C" mass="12.01078"/>
  <Type name="1705" class="O" element="O" mass="15.99943"/>
  <Type name="1706" class="CT" element="C" mass="12.01078"/>
  <Type name="1707" class="H1" element="H" mass="1.007947"/>
  <Type name="1708" class="CT" element="C" mass="12.01078"/>
  <Type name="1709" class="H1" element="H" mass="1.007947"/>
  <Type name="1710" class="OH" element="O" mass="15.99943"/>
  <Type name="1711" class="HO" element="H" mass="1.007947"/>
  <Type name="1712" class="OH" element="O" mass="15.99943"/>
  <Type name="1713" class="HO" element="H" mass="1.007947"/>
  <Type name="1714" class="P" element="P" mass="30.9737622"/>
  <Type name="1715" class="O2" element="O" mass="15.99943"/>
  <Type name="1716" class="O2" element="O" mass="15.99943"/>
  <Type name="1717" class="OS" element="O" mass="15.99943"/>
  <Type name="1718" class="CT" element="C" mass="12.01078"/>
  <Type name="1719" class="H1" element="H" mass="1.007947"/>
  <Type name="1720" class="CT" element="C" mass="12.01078"/>
  <Type name="1721" class="H1" element="H" mass="1.007947"/>
  <Type name="1722" class="OS" element="O" mass="15.99943"/>
  <Type name="1723" class="CT" element="C" mass="12.01078"/>
  <Type name="1724" class="H2" element="H" mass="1.007947"/>
  <Type name="1725" class="N*" element="N" mass="14.00672"/>
  <Type name="1726" class="CK" element="C" mass="12.01078"/>
  <Type name="1727" class="H5" element="H" mass="1.007947"/>
  <Type name="1728" class="NB" element="N" mass="14.00672"/>
  <Type name="1729" class="CB" element="C" mass="12.01078"/>
  <Type name="1730" class="C" element="C" mass="12.01078"/>
  <Type name="1731" class="O" element="O" mass="15.99943"/>
  <Type name="1732" class="NA" element="N" mass="14.00672"/>
  <Type name="1733" class="H" element="H" mass="1.007947"/>
  <Type name="1734" class="CA" element="C" mass="12.01078"/>
  <Type name="1735" class="N2" element="N" mass="14.00672"/>
  <Type name="1736" class="H" element="H" mass="1.007947"/>
  <Type name="1737" class="NC" element="N" mass="14.00672"/>
  <Type name="1738" class="CB" element="C" mass="12.01078"/>
  <Type name="1739" class="CT" element="C" mass="12.01078"/>
  <Type name="1740" class="H1" element="H" mass="1.007947"/>
  <Type name="1741" class="CT" element="C" mass="12.01078"/>
  <Type name="1742" class="H1" element="H" mass="1.007947"/>
  <Type name="1743" class="OH" element="O" mass="15.99943"/>
  <Type name="1744" class="HO" element="H" mass="1.007947"/>
  <Type name="1745" class="OS" element="O" mass="15.99943"/>
  <Type name="1746" class="P" element="P" mass="30.9737622"/>
  <Type name="1747" class="O2" element="O" mass="15.99943"/>
  <Type name="1748" class="O2" element="O" mass="15.99943"/>
  <Type name="1749" class="OS" element="O" mass="15.99943"/>
  <Type name="1750" class="CT" element="C" mass="12.01078"/>
  <Type name="1751" class="H1" element="H" mass="1.007947"/>
  <Type name="1752" class="CT" element="C" mass="12.01078"/>
  <Type name="1753" class="H1" element="H" mass="1.007947"/>
  <Type name="1754" class="OS" element="O" mass="15.99943"/>
  <Type name="1755" class="CT" element="C" mass="12.01078"/>
  <Type name="1756" class="H2" element="H" mass="1.007947"/>
  <Type name="1757" class="N*" element="N" mass="14.00672"/>
  <Type name="1758" class="CK" element="C" mass="12.01078"/>
  <Type name="1759" class="H5" element="H" mass="1.007947"/>
  <Type name="1760" class="NB" element="N" mass="14.00672"/>
  <Type name="1761" class="CB" element="C" mass="12.01078"/>
  <Type name="1762" class="C" element="C" mass="12.01078"/>
  <Type name="1763" class="O" element="O" mass="15.99943"/>
  <Type name="1764" class="NA" element="N" mass="14.00672"/>
  <Type name="1765" class="H" element="H" mass="1.007947"/>
  <Type name="1766" class="CA" element="C" mass="12.01078"/>
  <Type name="1767" class="N2" element="N" mass="14.00672"/>
  <Type name="1768" class="H" element="H" mass="1.007947"/>
  <Type name="1769" class="NC" element="N" mass="14.00672"/>
  <Type name="1770" class="CB" element="C" mass="12.01078"/>
  <Type name="1771" class="CT" element="C" mass="12.01078"/>
  <Type name="1772" class="H1" element="H" mass="1.007947"/>
  <Type name="1773" class="CT" element="C" mass="12.01078"/>
  <Type name="1774" class="H1" element="H" mass="1.007947"/>
  <Type name="1775" class="OH" element="O" mass="15.99943"/>
  <Type name="1776" class="HO" element="H" mass="1.007947"/>
  <Type name="1777" class="OH" element="O" mass="15.99943"/>
  <Type name="1778" class="HO" element="H" mass="1.007947"/>
  <Type name="1779" class="HO" element="H" mass="1.007947"/>
  <Type name="1780" class="OH" element="O" mass="15.99943"/>
  <Type name="1781" class="CT" element="C" mass="12.01078"/>
  <Type name="1782" class="H1" element="H" mass="1.007947"/>
  <Type name="1783" class="CT" element="C" mass="12.01078"/>
  <Type name="1784" class="H1" element="H" mass="1.007947"/>
  <Type name="1785" class="OS" element="O" mass="15.99943"/>
  <Type name="1786" class="CT" element="C" mass="12.01078"/>
  <Type name="1787" class="H2" element="H" mass="1.007947"/>
  <Type name="1788" class="N*" element="N" mass="14.00672"/>
  <Type name="1789" class="CK" element="C" mass="12.01078"/>
  <Type name="1790" class="H5" element="H" mass="1.007947"/>
  <Type name="1791" class="NB" element="N" mass="14.00672"/>
  <Type name="1792" class="CB" element="C" mass="12.01078"/>
  <Type name="1793" class="C" element="C" mass="12.01078"/>
  <Type name="1794" class="O" element="O" mass="15.99943"/>
  <Type name="1795" class="NA" element="N" mass="14.00672"/>
  <Type name="1796" class="H" element="H" mass="1.007947"/>
  <Type name="1797" class="CA" element="C" mass="12.01078"/>
  <Type name="1798" class="N2" element="N" mass="14.00672"/>
  <Type name="1799" class="H" element="H" mass="1.007947"/>
  <Type name="1800" class="NC" element="N" mass="14.00672"/>
  <Type name="1801" class="CB" element="C" mass="12.01078"/>
  <Type name="1802" class="CT" element="C" mass="12.01078"/>
  <Type name="1803" class="H1" element="H" mass="1.007947"/>
  <Type name="1804" class="CT" element="C" mass="12.01078"/>
  <Type name="1805" class="H1" element="H" mass="1.007947"/>
  <Type name="1806" class="OH" element="O" mass="15.99943"/>
  <Type name="1807" class="HO" element="H" mass="1.007947"/>
  <Type name="1808" class="OS" element="O" mass="15.99943"/>
  <Type name="1809" class="HO" element="H" mass="1.007947"/>
  <Type name="1810" class="OH" element="O" mass="15.99943"/>
  <Type name="1811" class="CT" element="C" mass="12.01078"/>
  <Type name="1812" class="H1" element="H" mass="1.007947"/>
  <Type name="1813" class="CT" element="C" mass="12.01078"/>
  <Type name="1814" class="H1" element="H" mass="1.007947"/>
  <Type name="1815" class="OS" element="O" mass="15.99943"/>
  <Type name="1816" class="CT" element="C" mass="12.01078"/>
  <Type name="1817" class="H2" element="H" mass="1.007947"/>
  <Type name="1818" class="N*" element="N" mass="14.00672"/>
  <Type name="1819" class="CK" element="C" mass="12.01078"/>
  <Type name="1820" class="H5" element="H" mass="1.007947"/>
  <Type name="1821" class="NB" element="N" mass="14.00672"/>
  <Type name="1822" class="CB" element="C" mass="12.01078"/>
  <Type name="1823" class="C" element="C" mass="12.01078"/>
  <Type name="1824" class="O" element="O" mass="15.99943"/>
  <Type name="1825" class="NA" element="N" mass="14.00672"/>
  <Type name="1826" class="H" element="H" mass="1.007947"/>
  <Type name="1827" class="CA" element="C" mass="12.01078"/>
  <Type name="1828" class="N2" element="N" mass="14.00672"/>
  <Type name="1829" class="H" element="H" mass="1.007947"/>
  <Type name="1830" class="NC" element="N" mass="14.00672"/>
  <Type name="1831" class="CB" element="C" mass="12.01078"/>
  <Type name="1832" class="CT" element="C" mass="12.01078"/>
  <Type name="1833" class="H1" element="H" mass="1.007947"/>
  <Type name="1834" class="CT" element="C" mass="12.01078"/>
  <Type name="1835" class="H1" element="H" mass="1.007947"/>
  <Type name="1836" class="OH" element="O" mass="15.99943"/>
  <Type name="1837" class="HO" element="H" mass="1.007947"/>
  <Type name="1838" class="OH" element="O" mass="15.99943"/>
  <Type name="1839" class="HO" element="H" mass="1.007947"/>
  <Type name="1840" class="P" element="P" mass="30.9737622"/>
  <Type name="1841" class="O2" element="O" mass="15.99943"/>
  <Type name="1842" class="O2" element="O" mass="15.99943"/>
  <Type name="1843" class="OS" element="O" mass="15.99943"/>
  <Type name="1844" class="CT" element="C" mass="12.01078"/>
  <Type name="1845" class="H1" element="H" mass="1.007947"/>
  <Type name="1846" class="CT" element="C" mass="12.01078"/>
  <Type name="1847" class="H1" element="H" mass="1.007947"/>
  <Type name="1848" class="OS" element="O" mass="15.99943"/>
  <Type name="1849" class="CT" element="C" mass="12.01078"/>
  <Type name="1850" class="H2" element="H" mass="1.007947"/>
  <Type name="1851" class="N*" element="N" mass="14.00672"/>
  <Type name="1852" class="CM" element="C" mass="12.01078"/>
  <Type name="1853" class="H4" element="H" mass="1.007947"/>
  <Type name="1854" class="CM" element="C" mass="12.01078"/>
  <Type name="1855" class="HA" element="H" mass="1.007947"/>
  <Type name="1856" class="C" element="C" mass="12.01078"/>
  <Type name="1857" class="O" element="O" mass="15.99943"/>
  <Type name="1858" class="NA" element="N" mass="14.00672"/>
  <Type name="1859" class="H" element="H" mass="1.007947"/>
  <Type name="1860" class="C" element="C" mass="12.01078"/>
  <Type name="1861" class="O" element="O" mass="15.99943"/>
  <Type name="1862" class="CT" element="C" mass="12.01078"/>
  <Type name="1863" class="H1" element="H" mass="1.007947"/>
  <Type name="1864" class="CT" element="C" mass="12.01078"/>
  <Type name="1865" class="H1" element="H" mass="1.007947"/>
  <Type name="1866" class="OH" element="O" mass="15.99943"/>
  <Type name="1867" class="HO" element="H" mass="1.007947"/>
  <Type name="1868" class="OS" element="O" mass="15.99943"/>
  <Type name="1869" class="P" element="P" mass="30.9737622"/>
  <Type name="1870" class="O2" element="O" mass="15.99943"/>
  <Type name="1871" class="O2" element="O" mass="15.99943"/>
  <Type name="1872" class="OS" element="O" mass="15.99943"/>
  <Type name="1873" class="CT" element="C" mass="12.01078"/>
  <Type name="1874" class="H1" element="H" mass="1.007947"/>
  <Type name="1875" class="CT" element="C" mass="12.01078"/>
  <Type name="1876" class="H1" element="H" mass="1.007947"/>
  <Type name="1877" class="OS" element="O" mass="15.99943"/>
  <Type name="1878" class="CT" element="C" mass="12.01078"/>
  <Type name="1879" class="H2" element="H" mass="1.007947"/>
  <Type name="1880" class="N*" element="N" mass="14.00672"/>
  <Type name="1881" class="CM" element="C" mass="12.01078"/>
  <Type name="1882" class="H4" element="H" mass="1.007947"/>
  <Type name="1883" class="CM" element="C" mass="12.01078"/>
  <Type name="1884" class="HA" element="H" mass="1.007947"/>
  <Type name="1885" class="C" element="C" mass="12.01078"/>
  <Type name="1886" class="O" element="O" mass="15.99943"/>
  <Type name="1887" class="NA" element="N" mass="14.00672"/>
  <Type name="1888" class="H" element="H" mass="1.007947"/>
  <Type name="1889" class="C" element="C" mass="12.01078"/>
  <Type name="1890" class="O" element="O" mass="15.99943"/>
  <Type name="1891" class="CT" element="C" mass="12.01078"/>
  <Type name="1892" class="H1" element="H" mass="1.007947"/>
  <Type name="1893" class="CT" element="C" mass="12.01078"/>
  <Type name="1894" class="H1" element="H" mass="1.007947"/>
  <Type name="1895" class="OH" element="O" mass="15.99943"/>
  <Type name="1896" class="HO" element="H" mass="1.007947"/>
  <Type name="1897" class="OH" element="O" mass="15.99943"/>
  <Type name="1898" class="HO" element="H" mass="1.007947"/>
  <Type name="1899" class="HO" element="H" mass="1.007947"/>
  <Type name="1900" class="OH" element="O" mass="15.99943"/>
  <Type name="1901" class="CT" element="C" mass="12.01078"/>
  <Type name="1902" class="H1" element="H" mass="1.007947"/>
  <Type name="1903" class="CT" element="C" mass="12.01078"/>
  <Type name="1904" class="H1" element="H" mass="1.007947"/>
  <Type name="1905" class="OS" element="O" mass="15.99943"/>
  <Type name="1906" class="CT" element="C" mass="12.01078"/>
  <Type name="1907" class="H2" element="H" mass="1.007947"/>
  <Type name="1908" class="N*" element="N" mass="14.00672"/>
  <Type name="1909" class="CM" element="C" mass="12.01078"/>
  <Type name="1910" class="H4" element="H" mass="1.007947"/>
  <Type name="1911" class="CM" element="C" mass="12.01078"/>
  <Type name="1912" class="HA" element="H" mass="1.007947"/>
  <Type name="1913" class="C" element="C" mass="12.01078"/>
  <Type name="1914" class="O" element="O" mass="15.99943"/>
  <Type name="1915" class="NA" element="N" mass="14.00672"/>
  <Type name="1916" class="H" element="H" mass="1.007947"/>
  <Type name="1917" class="C" element="C" mass="12.01078"/>
  <Type name="1918" class="O" element="O" mass="15.99943"/>
  <Type name="1919" class="CT" element="C" mass="12.01078"/>
  <Type name="1920" class="H1" element="H" mass="1.007947"/>
  <Type name="1921" class="CT" element="C" mass="12.01078"/>
  <Type name="1922" class="H1" element="H" mass="1.007947"/>
  <Type name="1923" class="OH" element="O" mass="15.99943"/>
  <Type name="1924" class="HO" element="H" mass="1.007947"/>
  <Type name="1925" class="OS" element="O" mass="15.99943"/>
  <Type name="1926" class="HO" element="H" mass="1.007947"/>
  <Type name="1927" class="OH" element="O" mass="15.99943"/>
  <Type name="1928" class="CT" element="C" mass="12.01078"/>
  <Type name="1929" class="H1" element="H" mass="1.007947"/>
  <Type name="1930" class="CT" element="C" mass="12.01078"/>
  <Type name="1931" class="H1" element="H" mass="1.007947"/>
  <Type name="1932" class="OS" element="O" mass="15.99943"/>
  <Type name="1933" class="CT" element="C" mass="12.01078"/>
  <Type name="1934" class="H2" element="H" mass="1.007947"/>
  <Type name="1935" class="N*" element="N" mass="14.00672"/>
  <Type name="1936" class="CM" element="C" mass="12.01078"/>
  <Type name="1937" class="H4" element="H" mass="1.007947"/>
  <Type name="1938" class="CM" element="C" mass="12.01078"/>
  <Type name="1939" class="HA" element="H" mass="1.007947"/>
  <Type name="1940" class="C" element="C" mass="12.01078"/>
  <Type name="1941" class="O" element="O" mass="15.99943"/>
  <Type name="1942" class="NA" element="N" mass="14.00672"/>
  <Type name="1943" class="H" element="H" mass="1.007947"/>
  <Type name="1944" class="C" element="C" mass="12.01078"/>
  <Type name="1945" class="O" element="O" mass="15.99943"/>
  <Type name="1946" class="CT" element="C" mass="12.01078"/>
  <Type name="1947" class="H1" element="H" mass="1.007947"/>
  <Type name="1948" class="CT" element="C" mass="12.01078"/>
  <Type name="1949" class="H1" element="H" mass="1.007947"/>
  <Type name="1950" class="OH" element="O" mass="15.99943"/>
  <Type name="1951" class="HO" element="H" mass="1.007947"/>
  <Type name="1952" class="OH" element="O" mass="15.99943"/>
  <Type name="1953" class="HO" element="H" mass="1.007947"/>
  <Type name="1954" class="IM" element="Cl" mass="35.4532"/>
  <Type name="1955" class="Cs" element="Cs" mass="132.90545192"/>
  <Type name="1956" class="K" element="K" mass="39.09831"/>
  <Type name="1957" class="Li" element="Li" mass="6.9412"/>
  <Type name="1958" class="MG" element="Mg" mass="24.30506"/>
  <Type name="1959" class="IP" element="Na" mass="22.989769282"/>
  <Type name="1960" class="Rb" element="Rb" mass="85.46783"/>
 </AtomTypes>
 <Residues>
  <Residue name="ACE">
   <Atom name="HH31" type="710"/>
   <Atom name="CH3" type="711"/>
   <Atom name="HH32" type="710"/>
   <Atom name="HH33" type="710"/>
   <Atom name="C" type="712"/>
   <Atom name="O" type="713"/>
   <Bond from="0" to="1"/>
   <Bond from="1" to="2"/>
   <Bond from="1" to="3"/>
   <Bond from="1" to="4"/>
   <Bond from="4" to="5"/>
   <ExternalBond from="4"/>
  </Residue>
  <Residue name="ALA">
   <Atom name="N" type="0"/>
   <Atom name="H" type="1"/>
   <Atom name="CA" type="2"/>
   <Atom name="HA" type="3"/>
   <Atom name="CB" type="4"/>
   <Atom name="HB1" type="5"/>
   <Atom name="HB2" type="5"/>
   <Atom name="HB3" type="5"/>
   <Atom name="C" type="6"/>
   <Atom name="O" type="7"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="8"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="8" to="9"/>
   <ExternalBond from="0"/>
   <ExternalBond from="8"/>
  </Residue>
  <Residue name="ARG">
   <Atom name="N" type="8"/>
   <Atom name="H" type="9"/>
   <Atom name="CA" type="10"/>
   <Atom name="HA" type="11"/>
   <Atom name="CB" type="12"/>
   <Atom name="HB2" type="13"/>
   <Atom name="HB3" type="13"/>
   <Atom name="CG" type="14"/>
   <Atom name="HG2" type="15"/>
   <Atom name="HG3" type="15"/>
   <Atom name="CD" type="16"/>
   <Atom name="HD2" type="17"/>
   <Atom name="HD3" type="17"/>
   <Atom name="NE" type="18"/>
   <Atom name="HE" type="19"/>
   <Atom name="CZ" type="20"/>
   <Atom name="NH1" type="21"/>
   <Atom name="HH11" type="22"/>
   <Atom name="HH12" type="22"/>
   <Atom name="NH2" type="23"/>
   <Atom name="HH21" type="24"/>
   <Atom name="HH22" type="24"/>
   <Atom name="C" type="25"/>
   <Atom name="O" type="26"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="22"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="7" to="9"/>
   <Bond from="7" to="10"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="12"/>
   <Bond from="10" to="13"/>
   <Bond from="13" to="14"/>
   <Bond from="13" to="15"/>
   <Bond from="15" to="16"/>
   <Bond from="15" to="19"/>
   <Bond from="16" to="17"/>
   <Bond from="16" to="18"/>
   <Bond from="19" to="20"/>
   <Bond from="19" to="21"/>
   <Bond from="22" to="23"/>
   <ExternalBond from="0"/>
   <ExternalBond from="22"/>
  </Residue>
  <Residue name="ASH">
   <Atom name="N" type="27"/>
   <Atom name="H" type="28"/>
   <Atom name="CA" type="29"/>
   <Atom name="HA" type="30"/>
   <Atom name="CB" type="31"/>
   <Atom name="HB2" type="32"/>
   <Atom name="HB3" type="32"/>
   <Atom name="CG" type="33"/>
   <Atom name="OD1" type="34"/>
   <Atom name="OD2" type="35"/>
   <Atom name="HD2" type="36"/>
   <Atom name="C" type="37"/>
   <Atom name="O" type="38"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="11"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="7" to="9"/>
   <Bond from="9" to="10"/>
   <Bond from="11" to="12"/>
   <ExternalBond from="0"/>
   <ExternalBond from="11"/>
  </Residue>
  <Residue name="ASN">
   <Atom name="N" type="39"/>
   <Atom name="H" type="40"/>
   <Atom name="CA" type="41"/>
   <Atom name="HA" type="42"/>
   <Atom name="CB" type="43"/>
   <Atom name="HB2" type="44"/>
   <Atom name="HB3" type="44"/>
   <Atom name="CG" type="45"/>
   <Atom name="OD1" type="46"/>
   <Atom name="ND2" type="47"/>
   <Atom name="HD21" type="48"/>
   <Atom name="HD22" type="48"/>
   <Atom name="C" type="49"/>
   <Atom name="O" type="50"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="12"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="7" to="9"/>
   <Bond from="9" to="10"/>
   <Bond from="9" to="11"/>
   <Bond from="12" to="13"/>
   <ExternalBond from="0"/>
   <ExternalBond from="12"/>
  </Residue>
  <Residue name="ASP">
   <Atom name="N" type="51"/>
   <Atom name="H" type="52"/>
   <Atom name="CA" type="53"/>
   <Atom name="HA" type="54"/>
   <Atom name="CB" type="55"/>
   <Atom name="HB2" type="56"/>
   <Atom name="HB3" type="56"/>
   <Atom name="CG" type="57"/>
   <Atom name="OD1" type="58"/>
   <Atom name="OD2" type="59"/>
   <Atom name="C" type="60"/>
   <Atom name="O" type="61"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="10"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="7" to="9"/>
   <Bond from="10" to="11"/>
   <ExternalBond from="0"/>
   <ExternalBond from="10"/>
  </Residue>
  <Residue name="CALA">
   <Atom name="N" type="366"/>
   <Atom name="H" type="367"/>
   <Atom name="CA" type="368"/>
   <Atom name="HA" type="369"/>
   <Atom name="CB" type="370"/>
   <Atom name="HB1" type="371"/>
   <Atom name="HB2" type="371"/>
   <Atom name="HB3" type="371"/>
   <Atom name="C" type="372"/>
   <Atom name="O" type="373"/>
   <Atom name="OXT" type="374"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="8"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="8" to="9"/>
   <Bond from="8" to="10"/>
   <ExternalBond from="0"/>
  </Residue>
  <Residue name="CARG">
   <Atom name="N" type="375"/>
   <Atom name="H" type="376"/>
   <Atom name="CA" type="377"/>
   <Atom name="HA" type="378"/>
   <Atom name="CB" type="379"/>
   <Atom name="HB2" type="380"/>
   <Atom name="HB3" type="380"/>
   <Atom name="CG" type="381"/>
   <Atom name="HG2" type="382"/>
   <Atom name="HG3" type="382"/>
   <Atom name="CD" type="383"/>
   <Atom name="HD2" type="384"/>
   <Atom name="HD3" type="384"/>
   <Atom name="NE" type="385"/>
   <Atom name="HE" type="386"/>
   <Atom name="CZ" type="387"/>
   <Atom name="NH1" type="388"/>
   <Atom name="HH11" type="389"/>
   <Atom name="HH12" type="389"/>
   <Atom name="NH2" type="390"/>
   <Atom name="HH21" type="391"/>
   <Atom name="HH22" type="391"/>
   <Atom name="C" type="392"/>
   <Atom name="O" type="393"/>
   <Atom name="OXT" type="394"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="22"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="7" to="9"/>
   <Bond from="7" to="10"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="12"/>
   <Bond from="10" to="13"/>
   <Bond from="13" to="14"/>
   <Bond from="13" to="15"/>
   <Bond from="15" to="16"/>
   <Bond from="15" to="19"/>
   <Bond from="16" to="17"/>
   <Bond from="16" to="18"/>
   <Bond from="19" to="20"/>
   <Bond from="19" to="21"/>
   <Bond from="22" to="23"/>
   <Bond from="22" to="24"/>
   <ExternalBond from="0"/>
  </Residue>
  <Residue name="CASN">
   <Atom name="N" type="395"/>
   <Atom name="H" type="396"/>
   <Atom name="CA" type="397"/>
   <Atom name="HA" type="398"/>
   <Atom name="CB" type="399"/>
   <Atom name="HB2" type="400"/>
   <Atom name="HB3" type="400"/>
   <Atom name="CG" type="401"/>
   <Atom name="OD1" type="402"/>
   <Atom name="ND2" type="403"/>
   <Atom name="HD21" type="404"/>
   <Atom name="HD22" type="404"/>
   <Atom name="C" type="405"/>
   <Atom name="O" type="406"/>
   <Atom name="OXT" type="407"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="12"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="7" to="9"/>
   <Bond from="9" to="10"/>
   <Bond from="9" to="11"/>
   <Bond from="12" to="13"/>
   <Bond from="12" to="14"/>
   <ExternalBond from="0"/>
  </Residue>
  <Residue name="CASP">
   <Atom name="N" type="408"/>
   <Atom name="H" type="409"/>
   <Atom name="CA" type="410"/>
   <Atom name="HA" type="411"/>
   <Atom name="CB" type="412"/>
   <Atom name="HB2" type="413"/>
   <Atom name="HB3" type="413"/>
   <Atom name="CG" type="414"/>
   <Atom name="OD1" type="415"/>
   <Atom name="OD2" type="416"/>
   <Atom name="C" type="417"/>
   <Atom name="O" type="418"/>
   <Atom name="OXT" type="419"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="10"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="7" to="9"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="12"/>
   <ExternalBond from="0"/>
  </Residue>
  <Residue name="CCYS">
   <Atom name="N" type="420"/>
   <Atom name="H" type="421"/>
   <Atom name="CA" type="422"/>
   <Atom name="HA" type="423"/>
   <Atom name="CB" type="424"/>
   <Atom name="HB2" type="425"/>
   <Atom name="HB3" type="425"/>
   <Atom name="SG" type="426"/>
   <Atom name="HG" type="427"/>
   <Atom name="C" type="428"/>
   <Atom name="O" type="429"/>
   <Atom name="OXT" type="430"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="9"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="9" to="10"/>
   <Bond from="9" to="11"/>
   <ExternalBond from="0"/>
  </Residue>
  <Residue name="CCYX">
   <Atom name="N" type="431"/>
   <Atom name="H" type="432"/>
   <Atom name="CA" type="433"/>
   <Atom name="HA" type="434"/>
   <Atom name="CB" type="435"/>
   <Atom name="HB2" type="436"/>
   <Atom name="HB3" type="436"/>
   <Atom name="SG" type="437"/>
   <Atom name="C" type="438"/>
   <Atom name="O" type="439"/>
   <Atom name="OXT" type="440"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="8"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="8" to="9"/>
   <Bond from="8" to="10"/>
   <ExternalBond from="0"/>
   <ExternalBond from="7"/>
  </Residue>
  <Residue name="CGLN">
   <Atom name="N" type="441"/>
   <Atom name="H" type="442"/>
   <Atom name="CA" type="443"/>
   <Atom name="HA" type="444"/>
   <Atom name="CB" type="445"/>
   <Atom name="HB2" type="446"/>
   <Atom name="HB3" type="446"/>
   <Atom name="CG" type="447"/>
   <Atom name="HG2" type="448"/>
   <Atom name="HG3" type="448"/>
   <Atom name="CD" type="449"/>
   <Atom name="OE1" type="450"/>
   <Atom name="NE2" type="451"/>
   <Atom name="HE21" type="452"/>
   <Atom name="HE22" type="452"/>
   <Atom name="C" type="453"/>
   <Atom name="O" type="454"/>
   <Atom name="OXT" type="455"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="15"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="7" to="9"/>
   <Bond from="7" to="10"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="12"/>
   <Bond from="12" to="13"/>
   <Bond from="12" to="14"/>
   <Bond from="15" to="16"/>
   <Bond from="15" to="17"/>
   <ExternalBond from="0"/>
  </Residue>
  <Residue name="CGLU">
   <Atom name="N" type="456"/>
   <Atom name="H" type="457"/>
   <Atom name="CA" type="458"/>
   <Atom name="HA" type="459"/>
   <Atom name="CB" type="460"/>
   <Atom name="HB2" type="461"/>
   <Atom name="HB3" type="461"/>
   <Atom name="CG" type="462"/>
   <Atom name="HG2" type="463"/>
   <Atom name="HG3" type="463"/>
   <Atom name="CD" type="464"/>
   <Atom name="OE1" type="465"/>
   <Atom name="OE2" type="466"/>
   <Atom name="C" type="467"/>
   <Atom name="O" type="468"/>
   <Atom name="OXT" type="469"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="13"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="7" to="9"/>
   <Bond from="7" to="10"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="12"/>
   <Bond from="13" to="14"/>
   <Bond from="13" to="15"/>
   <ExternalBond from="0"/>
  </Residue>
  <Residue name="CGLY">
   <Atom name="N" type="470"/>
   <Atom name="H" type="471"/>
   <Atom name="CA" type="472"/>
   <Atom name="HA2" type="473"/>
   <Atom name="HA3" type="473"/>
   <Atom name="C" type="474"/>
   <Atom name="O" type="475"/>
   <Atom name="OXT" type="476"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="5"/>
   <Bond from="5" to="6"/>
   <Bond from="5" to="7"/>
   <ExternalBond from="0"/>
  </Residue>
  <Residue name="CHID">
   <Atom name="N" type="477"/>
   <Atom name="H" type="478"/>
   <Atom name="CA" type="479"/>
   <Atom name="HA" type="480"/>
   <Atom name="CB" type="481"/>
   <Atom name="HB2" type="482"/>
   <Atom name="HB3" type="482"/>
   <Atom name="CG" type="483"/>
   <Atom name="ND1" type="484"/>
   <Atom name="HD1" type="485"/>
   <Atom name="CE1" type="486"/>
   <Atom name="HE1" type="487"/>
   <Atom name="NE2" type="488"/>
   <Atom name="CD2" type="489"/>
   <Atom name="HD2" type="490"/>
   <Atom name="C" type="491"/>
   <Atom name="O" type="492"/>
   <Atom name="OXT" type="493"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="15"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="7" to="13"/>
   <Bond from="8" to="9"/>
   <Bond from="8" to="10"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="12"/>
   <Bond from="12" to="13"/>
   <Bond from="13" to="14"/>
   <Bond from="15" to="16"/>
   <Bond from="15" to="17"/>
   <ExternalBond from="0"/>
  </Residue>
  <Residue name="CHIE">
   <Atom name="N" type="494"/>
   <Atom name="H" type="495"/>
   <Atom name="CA" type="496"/>
   <Atom name="HA" type="497"/>
   <Atom name="CB" type="498"/>
   <Atom name="HB2" type="499"/>
   <Atom name="HB3" type="499"/>
   <Atom name="CG" type="500"/>
   <Atom name="ND1" type="501"/>
   <Atom name="CE1" type="502"/>
   <Atom name="HE1" type="503"/>
   <Atom name="NE2" type="504"/>
   <Atom name="HE2" type="505"/>
   <Atom name="CD2" type="506"/>
   <Atom name="HD2" type="507"/>
   <Atom name="C" type="508"/>
   <Atom name="O" type="509"/>
   <Atom name="OXT" type="510"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="15"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="7" to="13"/>
   <Bond from="8" to="9"/>
   <Bond from="9" to="10"/>
   <Bond from="9" to="11"/>
   <Bond from="11" to="12"/>
   <Bond from="11" to="13"/>
   <Bond from="13" to="14"/>
   <Bond from="15" to="16"/>
   <Bond from="15" to="17"/>
   <ExternalBond from="0"/>
  </Residue>
  <Residue name="CHIP">
   <Atom name="N" type="511"/>
   <Atom name="H" type="512"/>
   <Atom name="CA" type="513"/>
   <Atom name="HA" type="514"/>
   <Atom name="CB" type="515"/>
   <Atom name="HB2" type="516"/>
   <Atom name="HB3" type="516"/>
   <Atom name="CG" type="517"/>
   <Atom name="ND1" type="518"/>
   <Atom name="HD1" type="519"/>
   <Atom name="CE1" type="520"/>
   <Atom name="HE1" type="521"/>
   <Atom name="NE2" type="522"/>
   <Atom name="HE2" type="523"/>
   <Atom name="CD2" type="524"/>
   <Atom name="HD2" type="525"/>
   <Atom name="C" type="526"/>
   <Atom name="O" type="527"/>
   <Atom name="OXT" type="528"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="16"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="7" to="14"/>
   <Bond from="8" to="9"/>
   <Bond from="8" to="10"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="12"/>
   <Bond from="12" to="13"/>
   <Bond from="12" to="14"/>
   <Bond from="14" to="15"/>
   <Bond from="16" to="17"/>
   <Bond from="16" to="18"/>
   <ExternalBond from="0"/>
  </Residue>
  <Residue name="CILE">
   <Atom name="N" type="529"/>
   <Atom name="H" type="530"/>
   <Atom name="CA" type="531"/>
   <Atom name="HA" type="532"/>
   <Atom name="CB" type="533"/>
   <Atom name="HB" type="534"/>
   <Atom name="CG2" type="535"/>
   <Atom name="HG21" type="536"/>
   <Atom name="HG22" type="536"/>
   <Atom name="HG23" type="536"/>
   <Atom name="CG1" type="537"/>
   <Atom name="HG12" type="538"/>
   <Atom name="HG13" type="538"/>
   <Atom name="CD1" type="539"/>
   <Atom name="HD11" type="540"/>
   <Atom name="HD12" type="540"/>
   <Atom name="HD13" type="540"/>
   <Atom name="C" type="541"/>
   <Atom name="O" type="542"/>
   <Atom name="OXT" type="543"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="17"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="10"/>
   <Bond from="6" to="7"/>
   <Bond from="6" to="8"/>
   <Bond from="6" to="9"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="12"/>
   <Bond from="10" to="13"/>
   <Bond from="13" to="14"/>
   <Bond from="13" to="15"/>
   <Bond from="13" to="16"/>
   <Bond from="17" to="18"/>
   <Bond from="17" to="19"/>
   <ExternalBond from="0"/>
  </Residue>
  <Residue name="CLEU">
   <Atom name="N" type="544"/>
   <Atom name="H" type="545"/>
   <Atom name="CA" type="546"/>
   <Atom name="HA" type="547"/>
   <Atom name="CB" type="548"/>
   <Atom name="HB2" type="549"/>
   <Atom name="HB3" type="549"/>
   <Atom name="CG" type="550"/>
   <Atom name="HG" type="551"/>
   <Atom name="CD1" type="552"/>
   <Atom name="HD11" type="553"/>
   <Atom name="HD12" type="553"/>
   <Atom name="HD13" type="553"/>
   <Atom name="CD2" type="554"/>
   <Atom name="HD21" type="555"/>
   <Atom name="HD22" type="555"/>
   <Atom name="HD23" type="555"/>
   <Atom name="C" type="556"/>
   <Atom name="O" type="557"/>
   <Atom name="OXT" type="558"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="17"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="7" to="9"/>
   <Bond from="7" to="13"/>
   <Bond from="9" to="10"/>
   <Bond from="9" to="11"/>
   <Bond from="9" to="12"/>
   <Bond from="13" to="14"/>
   <Bond from="13" to="15"/>
   <Bond from="13" to="16"/>
   <Bond from="17" to="18"/>
   <Bond from="17" to="19"/>
   <ExternalBond from="0"/>
  </Residue>
  <Residue name="CLYS">
   <Atom name="N" type="559"/>
   <Atom name="H" type="560"/>
   <Atom name="CA" type="561"/>
   <Atom name="HA" type="562"/>
   <Atom name="CB" type="563"/>
   <Atom name="HB2" type="564"/>
   <Atom name="HB3" type="564"/>
   <Atom name="CG" type="565"/>
   <Atom name="HG2" type="566"/>
   <Atom name="HG3" type="566"/>
   <Atom name="CD" type="567"/>
   <Atom name="HD2" type="568"/>
   <Atom name="HD3" type="568"/>
   <Atom name="CE" type="569"/>
   <Atom name="HE2" type="570"/>
   <Atom name="HE3" type="570"/>
   <Atom name="NZ" type="571"/>
   <Atom name="HZ1" type="572"/>
   <Atom name="HZ2" type="572"/>
   <Atom name="HZ3" type="572"/>
   <Atom name="C" type="573"/>
   <Atom name="O" type="574"/>
   <Atom name="OXT" type="575"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="20"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="7" to="9"/>
   <Bond from="7" to="10"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="12"/>
   <Bond from="10" to="13"/>
   <Bond from="13" to="14"/>
   <Bond from="13" to="15"/>
   <Bond from="13" to="16"/>
   <Bond from="16" to="17"/>
   <Bond from="16" to="18"/>
   <Bond from="16" to="19"/>
   <Bond from="20" to="21"/>
   <Bond from="20" to="22"/>
   <ExternalBond from="0"/>
  </Residue>
  <Residue name="CMET">
   <Atom name="N" type="576"/>
   <Atom name="H" type="577"/>
   <Atom name="CA" type="578"/>
   <Atom name="HA" type="579"/>
   <Atom name="CB" type="580"/>
   <Atom name="HB2" type="581"/>
   <Atom name="HB3" type="581"/>
   <Atom name="CG" type="582"/>
   <Atom name="HG2" type="583"/>
   <Atom name="HG3" type="583"/>
   <Atom name="SD" type="584"/>
   <Atom name="CE" type="585"/>
   <Atom name="HE1" type="586"/>
   <Atom name="HE2" type="586"/>
   <Atom name="HE3" type="586"/>
   <Atom name="C" type="587"/>
   <Atom name="O" type="588"/>
   <Atom name="OXT" type="589"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="15"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="7" to="9"/>
   <Bond from="7" to="10"/>
   <Bond from="10" to="11"/>
   <Bond from="11" to="12"/>
   <Bond from="11" to="13"/>
   <Bond from="11" to="14"/>
   <Bond from="15" to="16"/>
   <Bond from="15" to="17"/>
   <ExternalBond from="0"/>
  </Residue>
  <Residue name="CPHE">
   <Atom name="N" type="590"/>
   <Atom name="H" type="591"/>
   <Atom name="CA" type="592"/>
   <Atom name="HA" type="593"/>
   <Atom name="CB" type="594"/>
   <Atom name="HB2" type="595"/>
   <Atom name="HB3" type="595"/>
   <Atom name="CG" type="596"/>
   <Atom name="CD1" type="597"/>
   <Atom name="HD1" type="598"/>
   <Atom name="CE1" type="599"/>
   <Atom name="HE1" type="600"/>
   <Atom name="CZ" type="601"/>
   <Atom name="HZ" type="602"/>
   <Atom name="CE2" type="603"/>
   <Atom name="HE2" type="604"/>
   <Atom name="CD2" type="605"/>
   <Atom name="HD2" type="606"/>
   <Atom name="C" type="607"/>
   <Atom name="O" type="608"/>
   <Atom name="OXT" type="609"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="18"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="7" to="16"/>
   <Bond from="8" to="9"/>
   <Bond from="8" to="10"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="12"/>
   <Bond from="12" to="13"/>
   <Bond from="12" to="14"/>
   <Bond from="14" to="15"/>
   <Bond from="14" to="16"/>
   <Bond from="16" to="17"/>
   <Bond from="18" to="19"/>
   <Bond from="18" to="20"/>
   <ExternalBond from="0"/>
  </Residue>
  <Residue name="CPRO">
   <Atom name="N" type="610"/>
   <Atom name="CD" type="611"/>
   <Atom name="HD2" type="612"/>
   <Atom name="HD3" type="612"/>
   <Atom name="CG" type="613"/>
   <Atom name="HG2" type="614"/>
   <Atom name="HG3" type="614"/>
   <Atom name="CB" type="615"/>
   <Atom name="HB2" type="616"/>
   <Atom name="HB3" type="616"/>
   <Atom name="CA" type="617"/>
   <Atom name="HA" type="618"/>
   <Atom name="C" type="619"/>
   <Atom name="O" type="620"/>
   <Atom name="OXT" type="621"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="10"/>
   <Bond from="1" to="2"/>
   <Bond from="1" to="3"/>
   <Bond from="1" to="4"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="7" to="9"/>
   <Bond from="7" to="10"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="12"/>
   <Bond from="12" to="13"/>
   <Bond from="12" to="14"/>
   <ExternalBond from="0"/>
  </Residue>
  <Residue name="CSER">
   <Atom name="N" type="622"/>
   <Atom name="H" type="623"/>
   <Atom name="CA" type="624"/>
   <Atom name="HA" type="625"/>
   <Atom name="CB" type="626"/>
   <Atom name="HB2" type="627"/>
   <Atom name="HB3" type="627"/>
   <Atom name="OG" type="628"/>
   <Atom name="HG" type="629"/>
   <Atom name="C" type="630"/>
   <Atom name="O" type="631"/>
   <Atom name="OXT" type="632"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="9"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="9" to="10"/>
   <Bond from="9" to="11"/>
   <ExternalBond from="0"/>
  </Residue>
  <Residue name="CTHR">
   <Atom name="N" type="633"/>
   <Atom name="H" type="634"/>
   <Atom name="CA" type="635"/>
   <Atom name="HA" type="636"/>
   <Atom name="CB" type="637"/>
   <Atom name="HB" type="638"/>
   <Atom name="CG2" type="639"/>
   <Atom name="HG21" type="640"/>
   <Atom name="HG22" type="640"/>
   <Atom name="HG23" type="640"/>
   <Atom name="OG1" type="641"/>
   <Atom name="HG1" type="642"/>
   <Atom name="C" type="643"/>
   <Atom name="O" type="644"/>
   <Atom name="OXT" type="645"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="12"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="10"/>
   <Bond from="6" to="7"/>
   <Bond from="6" to="8"/>
   <Bond from="6" to="9"/>
   <Bond from="10" to="11"/>
   <Bond from="12" to="13"/>
   <Bond from="12" to="14"/>
   <ExternalBond from="0"/>
  </Residue>
  <Residue name="CTRP">
   <Atom name="N" type="646"/>
   <Atom name="H" type="647"/>
   <Atom name="CA" type="648"/>
   <Atom name="HA" type="649"/>
   <Atom name="CB" type="650"/>
   <Atom name="HB2" type="651"/>
   <Atom name="HB3" type="651"/>
   <Atom name="CG" type="652"/>
   <Atom name="CD1" type="653"/>
   <Atom name="HD1" type="654"/>
   <Atom name="NE1" type="655"/>
   <Atom name="HE1" type="656"/>
   <Atom name="CE2" type="657"/>
   <Atom name="CZ2" type="658"/>
   <Atom name="HZ2" type="659"/>
   <Atom name="CH2" type="660"/>
   <Atom name="HH2" type="661"/>
   <Atom name="CZ3" type="662"/>
   <Atom name="HZ3" type="663"/>
   <Atom name="CE3" type="664"/>
   <Atom name="HE3" type="665"/>
   <Atom name="CD2" type="666"/>
   <Atom name="C" type="667"/>
   <Atom name="O" type="668"/>
   <Atom name="OXT" type="669"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="22"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="7" to="21"/>
   <Bond from="8" to="9"/>
   <Bond from="8" to="10"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="12"/>
   <Bond from="12" to="13"/>
   <Bond from="12" to="21"/>
   <Bond from="13" to="14"/>
   <Bond from="13" to="15"/>
   <Bond from="15" to="16"/>
   <Bond from="15" to="17"/>
   <Bond from="17" to="18"/>
   <Bond from="17" to="19"/>
   <Bond from="19" to="20"/>
   <Bond from="19" to="21"/>
   <Bond from="22" to="23"/>
   <Bond from="22" to="24"/>
   <ExternalBond from="0"/>
  </Residue>
  <Residue name="CTYR">
   <Atom name="N" type="670"/>
   <Atom name="H" type="671"/>
   <Atom name="CA" type="672"/>
   <Atom name="HA" type="673"/>
   <Atom name="CB" type="674"/>
   <Atom name="HB2" type="675"/>
   <Atom name="HB3" type="675"/>
   <Atom name="CG" type="676"/>
   <Atom name="CD1" type="677"/>
   <Atom name="HD1" type="678"/>
   <Atom name="CE1" type="679"/>
   <Atom name="HE1" type="680"/>
   <Atom name="CZ" type="681"/>
   <Atom name="OH" type="682"/>
   <Atom name="HH" type="683"/>
   <Atom name="CE2" type="684"/>
   <Atom name="HE2" type="685"/>
   <Atom name="CD2" type="686"/>
   <Atom name="HD2" type="687"/>
   <Atom name="C" type="688"/>
   <Atom name="O" type="689"/>
   <Atom name="OXT" type="690"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="19"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="7" to="17"/>
   <Bond from="8" to="9"/>
   <Bond from="8" to="10"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="12"/>
   <Bond from="12" to="13"/>
   <Bond from="12" to="15"/>
   <Bond from="13" to="14"/>
   <Bond from="15" to="16"/>
   <Bond from="15" to="17"/>
   <Bond from="17" to="18"/>
   <Bond from="19" to="20"/>
   <Bond from="19" to="21"/>
   <ExternalBond from="0"/>
  </Residue>
  <Residue name="CVAL">
   <Atom name="N" type="691"/>
   <Atom name="H" type="692"/>
   <Atom name="CA" type="693"/>
   <Atom name="HA" type="694"/>
   <Atom name="CB" type="695"/>
   <Atom name="HB" type="696"/>
   <Atom name="CG1" type="697"/>
   <Atom name="HG11" type="698"/>
   <Atom name="HG12" type="698"/>
   <Atom name="HG13" type="698"/>
   <Atom name="CG2" type="699"/>
   <Atom name="HG21" type="700"/>
   <Atom name="HG22" type="700"/>
   <Atom name="HG23" type="700"/>
   <Atom name="C" type="701"/>
   <Atom name="O" type="702"/>
   <Atom name="OXT" type="703"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="14"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="10"/>
   <Bond from="6" to="7"/>
   <Bond from="6" to="8"/>
   <Bond from="6" to="9"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="12"/>
   <Bond from="10" to="13"/>
   <Bond from="14" to="15"/>
   <Bond from="14" to="16"/>
   <ExternalBond from="0"/>
  </Residue>
  <Residue name="CYM">
   <Atom name="N" type="62"/>
   <Atom name="H" type="63"/>
   <Atom name="CA" type="64"/>
   <Atom name="HA" type="65"/>
   <Atom name="CB" type="66"/>
   <Atom name="HB3" type="67"/>
   <Atom name="HB2" type="67"/>
   <Atom name="SG" type="68"/>
   <Atom name="C" type="69"/>
   <Atom name="O" type="70"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="8"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="8" to="9"/>
   <ExternalBond from="0"/>
   <ExternalBond from="8"/>
  </Residue>
  <Residue name="CYS">
   <Atom name="N" type="71"/>
   <Atom name="H" type="72"/>
   <Atom name="CA" type="73"/>
   <Atom name="HA" type="74"/>
   <Atom name="CB" type="75"/>
   <Atom name="HB2" type="76"/>
   <Atom name="HB3" type="76"/>
   <Atom name="SG" type="77"/>
   <Atom name="HG" type="78"/>
   <Atom name="C" type="79"/>
   <Atom name="O" type="80"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="9"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="9" to="10"/>
   <ExternalBond from="0"/>
   <ExternalBond from="9"/>
  </Residue>
  <Residue name="CYX">
   <Atom name="N" type="81"/>
   <Atom name="H" type="82"/>
   <Atom name="CA" type="83"/>
   <Atom name="HA" type="84"/>
   <Atom name="CB" type="85"/>
   <Atom name="HB2" type="86"/>
   <Atom name="HB3" type="86"/>
   <Atom name="SG" type="87"/>
   <Atom name="C" type="88"/>
   <Atom name="O" type="89"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="8"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="8" to="9"/>
   <ExternalBond from="0"/>
   <ExternalBond from="8"/>
   <ExternalBond from="7"/>
  </Residue>
  <Residue name="Cl-">
   <Atom name="Cl-" type="1954"/>
  </Residue>
  <Residue name="Cs+">
   <Atom name="Cs+" type="1955"/>
  </Residue>
  <Residue name="DA">
   <Atom name="P" type="1030"/>
   <Atom name="O1P" type="1031"/>
   <Atom name="O2P" type="1032"/>
   <Atom name="O5'" type="1033"/>
   <Atom name="C5'" type="1034"/>
   <Atom name="H5'1" type="1035"/>
   <Atom name="H5'2" type="1035"/>
   <Atom name="C4'" type="1036"/>
   <Atom name="H4'" type="1037"/>
   <Atom name="O4'" type="1038"/>
   <Atom name="C1'" type="1039"/>
   <Atom name="H1'" type="1040"/>
   <Atom name="N9" type="1041"/>
   <Atom name="C8" type="1042"/>
   <Atom name="H8" type="1043"/>
   <Atom name="N7" type="1044"/>
   <Atom name="C5" type="1045"/>
   <Atom name="C6" type="1046"/>
   <Atom name="N6" type="1047"/>
   <Atom name="H61" type="1048"/>
   <Atom name="H62" type="1048"/>
   <Atom name="N1" type="1049"/>
   <Atom name="C2" type="1050"/>
   <Atom name="H2" type="1051"/>
   <Atom name="N3" type="1052"/>
   <Atom name="C4" type="1053"/>
   <Atom name="C3'" type="1054"/>
   <Atom name="H3'" type="1055"/>
   <Atom name="C2'" type="1056"/>
   <Atom name="H2'1" type="1057"/>
   <Atom name="H2'2" type="1057"/>
   <Atom name="O3'" type="1058"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="0" to="3"/>
   <Bond from="3" to="4"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="7" to="9"/>
   <Bond from="7" to="26"/>
   <Bond from="9" to="10"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="12"/>
   <Bond from="10" to="28"/>
   <Bond from="12" to="13"/>
   <Bond from="12" to="25"/>
   <Bond from="13" to="14"/>
   <Bond from="13" to="15"/>
   <Bond from="15" to="16"/>
   <Bond from="16" to="17"/>
   <Bond from="16" to="25"/>
   <Bond from="17" to="18"/>
   <Bond from="17" to="21"/>
   <Bond from="18" to="19"/>
   <Bond from="18" to="20"/>
   <Bond from="21" to="22"/>
   <Bond from="22" to="23"/>
   <Bond from="22" to="24"/>
   <Bond from="24" to="25"/>
   <Bond from="26" to="27"/>
   <Bond from="26" to="28"/>
   <Bond from="26" to="31"/>
   <Bond from="28" to="29"/>
   <Bond from="28" to="30"/>
   <ExternalBond from="0"/>
   <ExternalBond from="31"/>
  </Residue>
  <Residue name="DA3">
   <Atom name="P" type="1059"/>
   <Atom name="O1P" type="1060"/>
   <Atom name="O2P" type="1061"/>
   <Atom name="O5'" type="1062"/>
   <Atom name="C5'" type="1063"/>
   <Atom name="H5'1" type="1064"/>
   <Atom name="H5'2" type="1064"/>
   <Atom name="C4'" type="1065"/>
   <Atom name="H4'" type="1066"/>
   <Atom name="O4'" type="1067"/>
   <Atom name="C1'" type="1068"/>
   <Atom name="H1'" type="1069"/>
   <Atom name="N9" type="1070"/>
   <Atom name="C8" type="1071"/>
   <Atom name="H8" type="1072"/>
   <Atom name="N7" type="1073"/>
   <Atom name="C5" type="1074"/>
   <Atom name="C6" type="1075"/>
   <Atom name="N6" type="1076"/>
   <Atom name="H61" type="1077"/>
   <Atom name="H62" type="1077"/>
   <Atom name="N1" type="1078"/>
   <Atom name="C2" type="1079"/>
   <Atom name="H2" type="1080"/>
   <Atom name="N3" type="1081"/>
   <Atom name="C4" type="1082"/>
   <Atom name="C3'" type="1083"/>
   <Atom name="H3'" type="1084"/>
   <Atom name="C2'" type="1085"/>
   <Atom name="H2'1" type="1086"/>
   <Atom name="H2'2" type="1086"/>
   <Atom name="O3'" type="1087"/>
   <Atom name="H3T" type="1088"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="0" to="3"/>
   <Bond from="3" to="4"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="7" to="9"/>
   <Bond from="7" to="26"/>
   <Bond from="9" to="10"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="12"/>
   <Bond from="10" to="28"/>
   <Bond from="12" to="13"/>
   <Bond from="12" to="25"/>
   <Bond from="13" to="14"/>
   <Bond from="13" to="15"/>
   <Bond from="15" to="16"/>
   <Bond from="16" to="17"/>
   <Bond from="16" to="25"/>
   <Bond from="17" to="18"/>
   <Bond from="17" to="21"/>
   <Bond from="18" to="19"/>
   <Bond from="18" to="20"/>
   <Bond from="21" to="22"/>
   <Bond from="22" to="23"/>
   <Bond from="22" to="24"/>
   <Bond from="24" to="25"/>
   <Bond from="26" to="27"/>
   <Bond from="26" to="28"/>
   <Bond from="26" to="31"/>
   <Bond from="28" to="29"/>
   <Bond from="28" to="30"/>
   <Bond from="31" to="32"/>
   <ExternalBond from="0"/>
  </Residue>
  <Residue name="DA5">
   <Atom name="H5T" type="1089"/>
   <Atom name="O5'" type="1090"/>
   <Atom name="C5'" type="1091"/>
   <Atom name="H5'1" type="1092"/>
   <Atom name="H5'2" type="1092"/>
   <Atom name="C4'" type="1093"/>
   <Atom name="H4'" type="1094"/>
   <Atom name="O4'" type="1095"/>
   <Atom name="C1'" type="1096"/>
   <Atom name="H1'" type="1097"/>
   <Atom name="N9" type="1098"/>
   <Atom name="C8" type="1099"/>
   <Atom name="H8" type="1100"/>
   <Atom name="N7" type="1101"/>
   <Atom name="C5" type="1102"/>
   <Atom name="C6" type="1103"/>
   <Atom name="N6" type="1104"/>
   <Atom name="H61" type="1105"/>
   <Atom name="H62" type="1105"/>
   <Atom name="N1" type="1106"/>
   <Atom name="C2" type="1107"/>
   <Atom name="H2" type="1108"/>
   <Atom name="N3" type="1109"/>
   <Atom name="C4" type="1110"/>
   <Atom name="C3'" type="1111"/>
   <Atom name="H3'" type="1112"/>
   <Atom name="C2'" type="1113"/>
   <Atom name="H2'1" type="1114"/>
   <Atom name="H2'2" type="1114"/>
   <Atom name="O3'" type="1115"/>
   <Bond from="0" to="1"/>
   <Bond from="1" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="5"/>
   <Bond from="5" to="6"/>
   <Bond from="5" to="7"/>
   <Bond from="5" to="24"/>
   <Bond from="7" to="8"/>
   <Bond from="8" to="9"/>
   <Bond from="8" to="10"/>
   <Bond from="8" to="26"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="23"/>
   <Bond from="11" to="12"/>
   <Bond from="11" to="13"/>
   <Bond from="13" to="14"/>
   <Bond from="14" to="15"/>
   <Bond from="14" to="23"/>
   <Bond from="15" to="16"/>
   <Bond from="15" to="19"/>
   <Bond from="16" to="17"/>
   <Bond from="16" to="18"/>
   <Bond from="19" to="20"/>
   <Bond from="20" to="21"/>
   <Bond from="20" to="22"/>
   <Bond from="22" to="23"/>
   <Bond from="24" to="25"/>
   <Bond from="24" to="26"/>
   <Bond from="24" to="29"/>
   <Bond from="26" to="27"/>
   <Bond from="26" to="28"/>
   <ExternalBond from="29"/>
  </Residue>
  <Residue name="DAN">
   <Atom name="H5T" type="1116"/>
   <Atom name="O5'" type="1117"/>
   <Atom name="C5'" type="1118"/>
   <Atom name="H5'1" type="1119"/>
   <Atom name="H5'2" type="1119"/>
   <Atom name="C4'" type="1120"/>
   <Atom name="H4'" type="1121"/>
   <Atom name="O4'" type="1122"/>
   <Atom name="C1'" type="1123"/>
   <Atom name="H1'" type="1124"/>
   <Atom name="N9" type="1125"/>
   <Atom name="C8" type="1126"/>
   <Atom name="H8" type="1127"/>
   <Atom name="N7" type="1128"/>
   <Atom name="C5" type="1129"/>
   <Atom name="C6" type="1130"/>
   <Atom name="N6" type="1131"/>
   <Atom name="H61" type="1132"/>
   <Atom name="H62" type="1132"/>
   <Atom name="N1" type="1133"/>
   <Atom name="C2" type="1134"/>
   <Atom name="H2" type="1135"/>
   <Atom name="N3" type="1136"/>
   <Atom name="C4" type="1137"/>
   <Atom name="C3'" type="1138"/>
   <Atom name="H3'" type="1139"/>
   <Atom name="C2'" type="1140"/>
   <Atom name="H2'1" type="1141"/>
   <Atom name="H2'2" type="1141"/>
   <Atom name="O3'" type="1142"/>
   <Atom name="H3T" type="1143"/>
   <Bond from="0" to="1"/>
   <Bond from="1" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="5"/>
   <Bond from="5" to="6"/>
   <Bond from="5" to="7"/>
   <Bond from="5" to="24"/>
   <Bond from="7" to="8"/>
   <Bond from="8" to="9"/>
   <Bond from="8" to="10"/>
   <Bond from="8" to="26"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="23"/>
   <Bond from="11" to="12"/>
   <Bond from="11" to="13"/>
   <Bond from="13" to="14"/>
   <Bond from="14" to="15"/>
   <Bond from="14" to="23"/>
   <Bond from="15" to="16"/>
   <Bond from="15" to="19"/>
   <Bond from="16" to="17"/>
   <Bond from="16" to="18"/>
   <Bond from="19" to="20"/>
   <Bond from="20" to="21"/>
   <Bond from="20" to="22"/>
   <Bond from="22" to="23"/>
   <Bond from="24" to="25"/>
   <Bond from="24" to="26"/>
   <Bond from="24" to="29"/>
   <Bond from="26" to="27"/>
   <Bond from="26" to="28"/>
   <Bond from="29" to="30"/>
  </Residue>
  <Residue name="DC">
   <Atom name="P" type="1144"/>
   <Atom name="O1P" type="1145"/>
   <Atom name="O2P" type="1146"/>
   <Atom name="O5'" type="1147"/>
   <Atom name="C5'" type="1148"/>
   <Atom name="H5'1" type="1149"/>
   <Atom name="H5'2" type="1149"/>
   <Atom name="C4'" type="1150"/>
   <Atom name="H4'" type="1151"/>
   <Atom name="O4'" type="1152"/>
   <Atom name="C1'" type="1153"/>
   <Atom name="H1'" type="1154"/>
   <Atom name="N1" type="1155"/>
   <Atom name="C6" type="1156"/>
   <Atom name="H6" type="1157"/>
   <Atom name="C5" type="1158"/>
   <Atom name="H5" type="1159"/>
   <Atom name="C4" type="1160"/>
   <Atom name="N4" type="1161"/>
   <Atom name="H41" type="1162"/>
   <Atom name="H42" type="1162"/>
   <Atom name="N3" type="1163"/>
   <Atom name="C2" type="1164"/>
   <Atom name="O2" type="1165"/>
   <Atom name="C3'" type="1166"/>
   <Atom name="H3'" type="1167"/>
   <Atom name="C2'" type="1168"/>
   <Atom name="H2'1" type="1169"/>
   <Atom name="H2'2" type="1169"/>
   <Atom name="O3'" type="1170"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="0" to="3"/>
   <Bond from="3" to="4"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="7" to="9"/>
   <Bond from="7" to="24"/>
   <Bond from="9" to="10"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="12"/>
   <Bond from="10" to="26"/>
   <Bond from="12" to="13"/>
   <Bond from="12" to="22"/>
   <Bond from="13" to="14"/>
   <Bond from="13" to="15"/>
   <Bond from="15" to="16"/>
   <Bond from="15" to="17"/>
   <Bond from="17" to="18"/>
   <Bond from="17" to="21"/>
   <Bond from="18" to="19"/>
   <Bond from="18" to="20"/>
   <Bond from="21" to="22"/>
   <Bond from="22" to="23"/>
   <Bond from="24" to="25"/>
   <Bond from="24" to="26"/>
   <Bond from="24" to="29"/>
   <Bond from="26" to="27"/>
   <Bond from="26" to="28"/>
   <ExternalBond from="0"/>
   <ExternalBond from="29"/>
  </Residue>
  <Residue name="DC3">
   <Atom name="P" type="1171"/>
   <Atom name="O1P" type="1172"/>
   <Atom name="O2P" type="1173"/>
   <Atom name="O5'" type="1174"/>
   <Atom name="C5'" type="1175"/>
   <Atom name="H5'1" type="1176"/>
   <Atom name="H5'2" type="1176"/>
   <Atom name="C4'" type="1177"/>
   <Atom name="H4'" type="1178"/>
   <Atom name="O4'" type="1179"/>
   <Atom name="C1'" type="1180"/>
   <Atom name="H1'" type="1181"/>
   <Atom name="N1" type="1182"/>
   <Atom name="C6" type="1183"/>
   <Atom name="H6" type="1184"/>
   <Atom name="C5" type="1185"/>
   <Atom name="H5" type="1186"/>
   <Atom name="C4" type="1187"/>
   <Atom name="N4" type="1188"/>
   <Atom name="H41" type="1189"/>
   <Atom name="H42" type="1189"/>
   <Atom name="N3" type="1190"/>
   <Atom name="C2" type="1191"/>
   <Atom name="O2" type="1192"/>
   <Atom name="C3'" type="1193"/>
   <Atom name="H3'" type="1194"/>
   <Atom name="C2'" type="1195"/>
   <Atom name="H2'1" type="1196"/>
   <Atom name="H2'2" type="1196"/>
   <Atom name="O3'" type="1197"/>
   <Atom name="H3T" type="1198"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="0" to="3"/>
   <Bond from="3" to="4"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="7" to="9"/>
   <Bond from="7" to="24"/>
   <Bond from="9" to="10"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="12"/>
   <Bond from="10" to="26"/>
   <Bond from="12" to="13"/>
   <Bond from="12" to="22"/>
   <Bond from="13" to="14"/>
   <Bond from="13" to="15"/>
   <Bond from="15" to="16"/>
   <Bond from="15" to="17"/>
   <Bond from="17" to="18"/>
   <Bond from="17" to="21"/>
   <Bond from="18" to="19"/>
   <Bond from="18" to="20"/>
   <Bond from="21" to="22"/>
   <Bond from="22" to="23"/>
   <Bond from="24" to="25"/>
   <Bond from="24" to="26"/>
   <Bond from="24" to="29"/>
   <Bond from="26" to="27"/>
   <Bond from="26" to="28"/>
   <Bond from="29" to="30"/>
   <ExternalBond from="0"/>
  </Residue>
  <Residue name="DC5">
   <Atom name="H5T" type="1199"/>
   <Atom name="O5'" type="1200"/>
   <Atom name="C5'" type="1201"/>
   <Atom name="H5'1" type="1202"/>
   <Atom name="H5'2" type="1202"/>
   <Atom name="C4'" type="1203"/>
   <Atom name="H4'" type="1204"/>
   <Atom name="O4'" type="1205"/>
   <Atom name="C1'" type="1206"/>
   <Atom name="H1'" type="1207"/>
   <Atom name="N1" type="1208"/>
   <Atom name="C6" type="1209"/>
   <Atom name="H6" type="1210"/>
   <Atom name="C5" type="1211"/>
   <Atom name="H5" type="1212"/>
   <Atom name="C4" type="1213"/>
   <Atom name="N4" type="1214"/>
   <Atom name="H41" type="1215"/>
   <Atom name="H42" type="1215"/>
   <Atom name="N3" type="1216"/>
   <Atom name="C2" type="1217"/>
   <Atom name="O2" type="1218"/>
   <Atom name="C3'" type="1219"/>
   <Atom name="H3'" type="1220"/>
   <Atom name="C2'" type="1221"/>
   <Atom name="H2'1" type="1222"/>
   <Atom name="H2'2" type="1222"/>
   <Atom name="O3'" type="1223"/>
   <Bond from="0" to="1"/>
   <Bond from="1" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="5"/>
   <Bond from="5" to="6"/>
   <Bond from="5" to="7"/>
   <Bond from="5" to="22"/>
   <Bond from="7" to="8"/>
   <Bond from="8" to="9"/>
   <Bond from="8" to="10"/>
   <Bond from="8" to="24"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="20"/>
   <Bond from="11" to="12"/>
   <Bond from="11" to="13"/>
   <Bond from="13" to="14"/>
   <Bond from="13" to="15"/>
   <Bond from="15" to="16"/>
   <Bond from="15" to="19"/>
   <Bond from="16" to="17"/>
   <Bond from="16" to="18"/>
   <Bond from="19" to="20"/>
   <Bond from="20" to="21"/>
   <Bond from="22" to="23"/>
   <Bond from="22" to="24"/>
   <Bond from="22" to="27"/>
   <Bond from="24" to="25"/>
   <Bond from="24" to="26"/>
   <ExternalBond from="27"/>
  </Residue>
  <Residue name="DCN">
   <Atom name="H5T" type="1224"/>
   <Atom name="O5'" type="1225"/>
   <Atom name="C5'" type="1226"/>
   <Atom name="H5'1" type="1227"/>
   <Atom name="H5'2" type="1227"/>
   <Atom name="C4'" type="1228"/>
   <Atom name="H4'" type="1229"/>
   <Atom name="O4'" type="1230"/>
   <Atom name="C1'" type="1231"/>
   <Atom name="H1'" type="1232"/>
   <Atom name="N1" type="1233"/>
   <Atom name="C6" type="1234"/>
   <Atom name="H6" type="1235"/>
   <Atom name="C5" type="1236"/>
   <Atom name="H5" type="1237"/>
   <Atom name="C4" type="1238"/>
   <Atom name="N4" type="1239"/>
   <Atom name="H41" type="1240"/>
   <Atom name="H42" type="1240"/>
   <Atom name="N3" type="1241"/>
   <Atom name="C2" type="1242"/>
   <Atom name="O2" type="1243"/>
   <Atom name="C3'" type="1244"/>
   <Atom name="H3'" type="1245"/>
   <Atom name="C2'" type="1246"/>
   <Atom name="H2'1" type="1247"/>
   <Atom name="H2'2" type="1247"/>
   <Atom name="O3'" type="1248"/>
   <Atom name="H3T" type="1249"/>
   <Bond from="0" to="1"/>
   <Bond from="1" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="5"/>
   <Bond from="5" to="6"/>
   <Bond from="5" to="7"/>
   <Bond from="5" to="22"/>
   <Bond from="7" to="8"/>
   <Bond from="8" to="9"/>
   <Bond from="8" to="10"/>
   <Bond from="8" to="24"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="20"/>
   <Bond from="11" to="12"/>
   <Bond from="11" to="13"/>
   <Bond from="13" to="14"/>
   <Bond from="13" to="15"/>
   <Bond from="15" to="16"/>
   <Bond from="15" to="19"/>
   <Bond from="16" to="17"/>
   <Bond from="16" to="18"/>
   <Bond from="19" to="20"/>
   <Bond from="20" to="21"/>
   <Bond from="22" to="23"/>
   <Bond from="22" to="24"/>
   <Bond from="22" to="27"/>
   <Bond from="24" to="25"/>
   <Bond from="24" to="26"/>
   <Bond from="27" to="28"/>
  </Residue>
  <Residue name="DG">
   <Atom name="P" type="1250"/>
   <Atom name="O1P" type="1251"/>
   <Atom name="O2P" type="1252"/>
   <Atom name="O5'" type="1253"/>
   <Atom name="C5'" type="1254"/>
   <Atom name="H5'1" type="1255"/>
   <Atom name="H5'2" type="1255"/>
   <Atom name="C4'" type="1256"/>
   <Atom name="H4'" type="1257"/>
   <Atom name="O4'" type="1258"/>
   <Atom name="C1'" type="1259"/>
   <Atom name="H1'" type="1260"/>
   <Atom name="N9" type="1261"/>
   <Atom name="C8" type="1262"/>
   <Atom name="H8" type="1263"/>
   <Atom name="N7" type="1264"/>
   <Atom name="C5" type="1265"/>
   <Atom name="C6" type="1266"/>
   <Atom name="O6" type="1267"/>
   <Atom name="N1" type="1268"/>
   <Atom name="H1" type="1269"/>
   <Atom name="C2" type="1270"/>
   <Atom name="N2" type="1271"/>
   <Atom name="H21" type="1272"/>
   <Atom name="H22" type="1272"/>
   <Atom name="N3" type="1273"/>
   <Atom name="C4" type="1274"/>
   <Atom name="C3'" type="1275"/>
   <Atom name="H3'" type="1276"/>
   <Atom name="C2'" type="1277"/>
   <Atom name="H2'1" type="1278"/>
   <Atom name="H2'2" type="1278"/>
   <Atom name="O3'" type="1279"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="0" to="3"/>
   <Bond from="3" to="4"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="7" to="9"/>
   <Bond from="7" to="27"/>
   <Bond from="9" to="10"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="12"/>
   <Bond from="10" to="29"/>
   <Bond from="12" to="13"/>
   <Bond from="12" to="26"/>
   <Bond from="13" to="14"/>
   <Bond from="13" to="15"/>
   <Bond from="15" to="16"/>
   <Bond from="16" to="17"/>
   <Bond from="16" to="26"/>
   <Bond from="17" to="18"/>
   <Bond from="17" to="19"/>
   <Bond from="19" to="20"/>
   <Bond from="19" to="21"/>
   <Bond from="21" to="22"/>
   <Bond from="21" to="25"/>
   <Bond from="22" to="23"/>
   <Bond from="22" to="24"/>
   <Bond from="25" to="26"/>
   <Bond from="27" to="28"/>
   <Bond from="27" to="29"/>
   <Bond from="27" to="32"/>
   <Bond from="29" to="30"/>
   <Bond from="29" to="31"/>
   <ExternalBond from="0"/>
   <ExternalBond from="32"/>
  </Residue>
  <Residue name="DG3">
   <Atom name="P" type="1280"/>
   <Atom name="O1P" type="1281"/>
   <Atom name="O2P" type="1282"/>
   <Atom name="O5'" type="1283"/>
   <Atom name="C5'" type="1284"/>
   <Atom name="H5'1" type="1285"/>
   <Atom name="H5'2" type="1285"/>
   <Atom name="C4'" type="1286"/>
   <Atom name="H4'" type="1287"/>
   <Atom name="O4'" type="1288"/>
   <Atom name="C1'" type="1289"/>
   <Atom name="H1'" type="1290"/>
   <Atom name="N9" type="1291"/>
   <Atom name="C8" type="1292"/>
   <Atom name="H8" type="1293"/>
   <Atom name="N7" type="1294"/>
   <Atom name="C5" type="1295"/>
   <Atom name="C6" type="1296"/>
   <Atom name="O6" type="1297"/>
   <Atom name="N1" type="1298"/>
   <Atom name="H1" type="1299"/>
   <Atom name="C2" type="1300"/>
   <Atom name="N2" type="1301"/>
   <Atom name="H21" type="1302"/>
   <Atom name="H22" type="1302"/>
   <Atom name="N3" type="1303"/>
   <Atom name="C4" type="1304"/>
   <Atom name="C3'" type="1305"/>
   <Atom name="H3'" type="1306"/>
   <Atom name="C2'" type="1307"/>
   <Atom name="H2'1" type="1308"/>
   <Atom name="H2'2" type="1308"/>
   <Atom name="O3'" type="1309"/>
   <Atom name="H3T" type="1310"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="0" to="3"/>
   <Bond from="3" to="4"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="7" to="9"/>
   <Bond from="7" to="27"/>
   <Bond from="9" to="10"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="12"/>
   <Bond from="10" to="29"/>
   <Bond from="12" to="13"/>
   <Bond from="12" to="26"/>
   <Bond from="13" to="14"/>
   <Bond from="13" to="15"/>
   <Bond from="15" to="16"/>
   <Bond from="16" to="17"/>
   <Bond from="16" to="26"/>
   <Bond from="17" to="18"/>
   <Bond from="17" to="19"/>
   <Bond from="19" to="20"/>
   <Bond from="19" to="21"/>
   <Bond from="21" to="22"/>
   <Bond from="21" to="25"/>
   <Bond from="22" to="23"/>
   <Bond from="22" to="24"/>
   <Bond from="25" to="26"/>
   <Bond from="27" to="28"/>
   <Bond from="27" to="29"/>
   <Bond from="27" to="32"/>
   <Bond from="29" to="30"/>
   <Bond from="29" to="31"/>
   <Bond from="32" to="33"/>
   <ExternalBond from="0"/>
  </Residue>
  <Residue name="DG5">
   <Atom name="H5T" type="1311"/>
   <Atom name="O5'" type="1312"/>
   <Atom name="C5'" type="1313"/>
   <Atom name="H5'1" type="1314"/>
   <Atom name="H5'2" type="1314"/>
   <Atom name="C4'" type="1315"/>
   <Atom name="H4'" type="1316"/>
   <Atom name="O4'" type="1317"/>
   <Atom name="C1'" type="1318"/>
   <Atom name="H1'" type="1319"/>
   <Atom name="N9" type="1320"/>
   <Atom name="C8" type="1321"/>
   <Atom name="H8" type="1322"/>
   <Atom name="N7" type="1323"/>
   <Atom name="C5" type="1324"/>
   <Atom name="C6" type="1325"/>
   <Atom name="O6" type="1326"/>
   <Atom name="N1" type="1327"/>
   <Atom name="H1" type="1328"/>
   <Atom name="C2" type="1329"/>
   <Atom name="N2" type="1330"/>
   <Atom name="H21" type="1331"/>
   <Atom name="H22" type="1331"/>
   <Atom name="N3" type="1332"/>
   <Atom name="C4" type="1333"/>
   <Atom name="C3'" type="1334"/>
   <Atom name="H3'" type="1335"/>
   <Atom name="C2'" type="1336"/>
   <Atom name="H2'1" type="1337"/>
   <Atom name="H2'2" type="1337"/>
   <Atom name="O3'" type="1338"/>
   <Bond from="0" to="1"/>
   <Bond from="1" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="5"/>
   <Bond from="5" to="6"/>
   <Bond from="5" to="7"/>
   <Bond from="5" to="25"/>
   <Bond from="7" to="8"/>
   <Bond from="8" to="9"/>
   <Bond from="8" to="10"/>
   <Bond from="8" to="27"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="24"/>
   <Bond from="11" to="12"/>
   <Bond from="11" to="13"/>
   <Bond from="13" to="14"/>
   <Bond from="14" to="15"/>
   <Bond from="14" to="24"/>
   <Bond from="15" to="16"/>
   <Bond from="15" to="17"/>
   <Bond from="17" to="18"/>
   <Bond from="17" to="19"/>
   <Bond from="19" to="20"/>
   <Bond from="19" to="23"/>
   <Bond from="20" to="21"/>
   <Bond from="20" to="22"/>
   <Bond from="23" to="24"/>
   <Bond from="25" to="26"/>
   <Bond from="25" to="27"/>
   <Bond from="25" to="30"/>
   <Bond from="27" to="28"/>
   <Bond from="27" to="29"/>
   <ExternalBond from="30"/>
  </Residue>
  <Residue name="DGN">
   <Atom name="H5T" type="1339"/>
   <Atom name="O5'" type="1340"/>
   <Atom name="C5'" type="1341"/>
   <Atom name="H5'1" type="1342"/>
   <Atom name="H5'2" type="1342"/>
   <Atom name="C4'" type="1343"/>
   <Atom name="H4'" type="1344"/>
   <Atom name="O4'" type="1345"/>
   <Atom name="C1'" type="1346"/>
   <Atom name="H1'" type="1347"/>
   <Atom name="N9" type="1348"/>
   <Atom name="C8" type="1349"/>
   <Atom name="H8" type="1350"/>
   <Atom name="N7" type="1351"/>
   <Atom name="C5" type="1352"/>
   <Atom name="C6" type="1353"/>
   <Atom name="O6" type="1354"/>
   <Atom name="N1" type="1355"/>
   <Atom name="H1" type="1356"/>
   <Atom name="C2" type="1357"/>
   <Atom name="N2" type="1358"/>
   <Atom name="H21" type="1359"/>
   <Atom name="H22" type="1359"/>
   <Atom name="N3" type="1360"/>
   <Atom name="C4" type="1361"/>
   <Atom name="C3'" type="1362"/>
   <Atom name="H3'" type="1363"/>
   <Atom name="C2'" type="1364"/>
   <Atom name="H2'1" type="1365"/>
   <Atom name="H2'2" type="1365"/>
   <Atom name="O3'" type="1366"/>
   <Atom name="H3T" type="1367"/>
   <Bond from="0" to="1"/>
   <Bond from="1" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="5"/>
   <Bond from="5" to="6"/>
   <Bond from="5" to="7"/>
   <Bond from="5" to="25"/>
   <Bond from="7" to="8"/>
   <Bond from="8" to="9"/>
   <Bond from="8" to="10"/>
   <Bond from="8" to="27"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="24"/>
   <Bond from="11" to="12"/>
   <Bond from="11" to="13"/>
   <Bond from="13" to="14"/>
   <Bond from="14" to="15"/>
   <Bond from="14" to="24"/>
   <Bond from="15" to="16"/>
   <Bond from="15" to="17"/>
   <Bond from="17" to="18"/>
   <Bond from="17" to="19"/>
   <Bond from="19" to="20"/>
   <Bond from="19" to="23"/>
   <Bond from="20" to="21"/>
   <Bond from="20" to="22"/>
   <Bond from="23" to="24"/>
   <Bond from="25" to="26"/>
   <Bond from="25" to="27"/>
   <Bond from="25" to="30"/>
   <Bond from="27" to="28"/>
   <Bond from="27" to="29"/>
   <Bond from="30" to="31"/>
  </Residue>
  <Residue name="DT">
   <Atom name="P" type="1368"/>
   <Atom name="O1P" type="1369"/>
   <Atom name="O2P" type="1370"/>
   <Atom name="O5'" type="1371"/>
   <Atom name="C5'" type="1372"/>
   <Atom name="H5'1" type="1373"/>
   <Atom name="H5'2" type="1373"/>
   <Atom name="C4'" type="1374"/>
   <Atom name="H4'" type="1375"/>
   <Atom name="O4'" type="1376"/>
   <Atom name="C1'" type="1377"/>
   <Atom name="H1'" type="1378"/>
   <Atom name="N1" type="1379"/>
   <Atom name="C6" type="1380"/>
   <Atom name="H6" type="1381"/>
   <Atom name="C5" type="1382"/>
   <Atom name="C7" type="1383"/>
   <Atom name="H71" type="1384"/>
   <Atom name="H72" type="1384"/>
   <Atom name="H73" type="1384"/>
   <Atom name="C4" type="1385"/>
   <Atom name="O4" type="1386"/>
   <Atom name="N3" type="1387"/>
   <Atom name="H3" type="1388"/>
   <Atom name="C2" type="1389"/>
   <Atom name="O2" type="1390"/>
   <Atom name="C3'" type="1391"/>
   <Atom name="H3'" type="1392"/>
   <Atom name="C2'" type="1393"/>
   <Atom name="H2'1" type="1394"/>
   <Atom name="H2'2" type="1394"/>
   <Atom name="O3'" type="1395"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="0" to="3"/>
   <Bond from="3" to="4"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="7" to="9"/>
   <Bond from="7" to="26"/>
   <Bond from="9" to="10"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="12"/>
   <Bond from="10" to="28"/>
   <Bond from="12" to="13"/>
   <Bond from="12" to="24"/>
   <Bond from="13" to="14"/>
   <Bond from="13" to="15"/>
   <Bond from="15" to="16"/>
   <Bond from="15" to="20"/>
   <Bond from="16" to="17"/>
   <Bond from="16" to="18"/>
   <Bond from="16" to="19"/>
   <Bond from="20" to="21"/>
   <Bond from="20" to="22"/>
   <Bond from="22" to="23"/>
   <Bond from="22" to="24"/>
   <Bond from="24" to="25"/>
   <Bond from="26" to="27"/>
   <Bond from="26" to="28"/>
   <Bond from="26" to="31"/>
   <Bond from="28" to="29"/>
   <Bond from="28" to="30"/>
   <ExternalBond from="0"/>
   <ExternalBond from="31"/>
  </Residue>
  <Residue name="DT3">
   <Atom name="P" type="1396"/>
   <Atom name="O1P" type="1397"/>
   <Atom name="O2P" type="1398"/>
   <Atom name="O5'" type="1399"/>
   <Atom name="C5'" type="1400"/>
   <Atom name="H5'1" type="1401"/>
   <Atom name="H5'2" type="1401"/>
   <Atom name="C4'" type="1402"/>
   <Atom name="H4'" type="1403"/>
   <Atom name="O4'" type="1404"/>
   <Atom name="C1'" type="1405"/>
   <Atom name="H1'" type="1406"/>
   <Atom name="N1" type="1407"/>
   <Atom name="C6" type="1408"/>
   <Atom name="H6" type="1409"/>
   <Atom name="C5" type="1410"/>
   <Atom name="C7" type="1411"/>
   <Atom name="H71" type="1412"/>
   <Atom name="H72" type="1412"/>
   <Atom name="H73" type="1412"/>
   <Atom name="C4" type="1413"/>
   <Atom name="O4" type="1414"/>
   <Atom name="N3" type="1415"/>
   <Atom name="H3" type="1416"/>
   <Atom name="C2" type="1417"/>
   <Atom name="O2" type="1418"/>
   <Atom name="C3'" type="1419"/>
   <Atom name="H3'" type="1420"/>
   <Atom name="C2'" type="1421"/>
   <Atom name="H2'1" type="1422"/>
   <Atom name="H2'2" type="1422"/>
   <Atom name="O3'" type="1423"/>
   <Atom name="H3T" type="1424"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="0" to="3"/>
   <Bond from="3" to="4"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="7" to="9"/>
   <Bond from="7" to="26"/>
   <Bond from="9" to="10"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="12"/>
   <Bond from="10" to="28"/>
   <Bond from="12" to="13"/>
   <Bond from="12" to="24"/>
   <Bond from="13" to="14"/>
   <Bond from="13" to="15"/>
   <Bond from="15" to="16"/>
   <Bond from="15" to="20"/>
   <Bond from="16" to="17"/>
   <Bond from="16" to="18"/>
   <Bond from="16" to="19"/>
   <Bond from="20" to="21"/>
   <Bond from="20" to="22"/>
   <Bond from="22" to="23"/>
   <Bond from="22" to="24"/>
   <Bond from="24" to="25"/>
   <Bond from="26" to="27"/>
   <Bond from="26" to="28"/>
   <Bond from="26" to="31"/>
   <Bond from="28" to="29"/>
   <Bond from="28" to="30"/>
   <Bond from="31" to="32"/>
   <ExternalBond from="0"/>
  </Residue>
  <Residue name="DT5">
   <Atom name="H5T" type="1425"/>
   <Atom name="O5'" type="1426"/>
   <Atom name="C5'" type="1427"/>
   <Atom name="H5'1" type="1428"/>
   <Atom name="H5'2" type="1428"/>
   <Atom name="C4'" type="1429"/>
   <Atom name="H4'" type="1430"/>
   <Atom name="O4'" type="1431"/>
   <Atom name="C1'" type="1432"/>
   <Atom name="H1'" type="1433"/>
   <Atom name="N1" type="1434"/>
   <Atom name="C6" type="1435"/>
   <Atom name="H6" type="1436"/>
   <Atom name="C5" type="1437"/>
   <Atom name="C7" type="1438"/>
   <Atom name="H71" type="1439"/>
   <Atom name="H72" type="1439"/>
   <Atom name="H73" type="1439"/>
   <Atom name="C4" type="1440"/>
   <Atom name="O4" type="1441"/>
   <Atom name="N3" type="1442"/>
   <Atom name="H3" type="1443"/>
   <Atom name="C2" type="1444"/>
   <Atom name="O2" type="1445"/>
   <Atom name="C3'" type="1446"/>
   <Atom name="H3'" type="1447"/>
   <Atom name="C2'" type="1448"/>
   <Atom name="H2'1" type="1449"/>
   <Atom name="H2'2" type="1449"/>
   <Atom name="O3'" type="1450"/>
   <Bond from="0" to="1"/>
   <Bond from="1" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="5"/>
   <Bond from="5" to="6"/>
   <Bond from="5" to="7"/>
   <Bond from="5" to="24"/>
   <Bond from="7" to="8"/>
   <Bond from="8" to="9"/>
   <Bond from="8" to="10"/>
   <Bond from="8" to="26"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="22"/>
   <Bond from="11" to="12"/>
   <Bond from="11" to="13"/>
   <Bond from="13" to="14"/>
   <Bond from="13" to="18"/>
   <Bond from="14" to="15"/>
   <Bond from="14" to="16"/>
   <Bond from="14" to="17"/>
   <Bond from="18" to="19"/>
   <Bond from="18" to="20"/>
   <Bond from="20" to="21"/>
   <Bond from="20" to="22"/>
   <Bond from="22" to="23"/>
   <Bond from="24" to="25"/>
   <Bond from="24" to="26"/>
   <Bond from="24" to="29"/>
   <Bond from="26" to="27"/>
   <Bond from="26" to="28"/>
   <ExternalBond from="29"/>
  </Residue>
  <Residue name="DTN">
   <Atom name="H5T" type="1451"/>
   <Atom name="O5'" type="1452"/>
   <Atom name="C5'" type="1453"/>
   <Atom name="H5'1" type="1454"/>
   <Atom name="H5'2" type="1454"/>
   <Atom name="C4'" type="1455"/>
   <Atom name="H4'" type="1456"/>
   <Atom name="O4'" type="1457"/>
   <Atom name="C1'" type="1458"/>
   <Atom name="H1'" type="1459"/>
   <Atom name="N1" type="1460"/>
   <Atom name="C6" type="1461"/>
   <Atom name="H6" type="1462"/>
   <Atom name="C5" type="1463"/>
   <Atom name="C7" type="1464"/>
   <Atom name="H71" type="1465"/>
   <Atom name="H72" type="1465"/>
   <Atom name="H73" type="1465"/>
   <Atom name="C4" type="1466"/>
   <Atom name="O4" type="1467"/>
   <Atom name="N3" type="1468"/>
   <Atom name="H3" type="1469"/>
   <Atom name="C2" type="1470"/>
   <Atom name="O2" type="1471"/>
   <Atom name="C3'" type="1472"/>
   <Atom name="H3'" type="1473"/>
   <Atom name="C2'" type="1474"/>
   <Atom name="H2'1" type="1475"/>
   <Atom name="H2'2" type="1475"/>
   <Atom name="O3'" type="1476"/>
   <Atom name="H3T" type="1477"/>
   <Bond from="0" to="1"/>
   <Bond from="1" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="5"/>
   <Bond from="5" to="6"/>
   <Bond from="5" to="7"/>
   <Bond from="5" to="24"/>
   <Bond from="7" to="8"/>
   <Bond from="8" to="9"/>
   <Bond from="8" to="10"/>
   <Bond from="8" to="26"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="22"/>
   <Bond from="11" to="12"/>
   <Bond from="11" to="13"/>
   <Bond from="13" to="14"/>
   <Bond from="13" to="18"/>
   <Bond from="14" to="15"/>
   <Bond from="14" to="16"/>
   <Bond from="14" to="17"/>
   <Bond from="18" to="19"/>
   <Bond from="18" to="20"/>
   <Bond from="20" to="21"/>
   <Bond from="20" to="22"/>
   <Bond from="22" to="23"/>
   <Bond from="24" to="25"/>
   <Bond from="24" to="26"/>
   <Bond from="24" to="29"/>
   <Bond from="26" to="27"/>
   <Bond from="26" to="28"/>
   <Bond from="29" to="30"/>
  </Residue>
  <Residue name="GLH">
   <Atom name="N" type="90"/>
   <Atom name="H" type="91"/>
   <Atom name="CA" type="92"/>
   <Atom name="HA" type="93"/>
   <Atom name="CB" type="94"/>
   <Atom name="HB2" type="95"/>
   <Atom name="HB3" type="95"/>
   <Atom name="CG" type="96"/>
   <Atom name="HG2" type="97"/>
   <Atom name="HG3" type="97"/>
   <Atom name="CD" type="98"/>
   <Atom name="OE1" type="99"/>
   <Atom name="OE2" type="100"/>
   <Atom name="HE2" type="101"/>
   <Atom name="C" type="102"/>
   <Atom name="O" type="103"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="14"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="7" to="9"/>
   <Bond from="7" to="10"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="12"/>
   <Bond from="12" to="13"/>
   <Bond from="14" to="15"/>
   <ExternalBond from="0"/>
   <ExternalBond from="14"/>
  </Residue>
  <Residue name="GLN">
   <Atom name="N" type="104"/>
   <Atom name="H" type="105"/>
   <Atom name="CA" type="106"/>
   <Atom name="HA" type="107"/>
   <Atom name="CB" type="108"/>
   <Atom name="HB2" type="109"/>
   <Atom name="HB3" type="109"/>
   <Atom name="CG" type="110"/>
   <Atom name="HG2" type="111"/>
   <Atom name="HG3" type="111"/>
   <Atom name="CD" type="112"/>
   <Atom name="OE1" type="113"/>
   <Atom name="NE2" type="114"/>
   <Atom name="HE21" type="115"/>
   <Atom name="HE22" type="115"/>
   <Atom name="C" type="116"/>
   <Atom name="O" type="117"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="15"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="7" to="9"/>
   <Bond from="7" to="10"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="12"/>
   <Bond from="12" to="13"/>
   <Bond from="12" to="14"/>
   <Bond from="15" to="16"/>
   <ExternalBond from="0"/>
   <ExternalBond from="15"/>
  </Residue>
  <Residue name="GLU">
   <Atom name="N" type="118"/>
   <Atom name="H" type="119"/>
   <Atom name="CA" type="120"/>
   <Atom name="HA" type="121"/>
   <Atom name="CB" type="122"/>
   <Atom name="HB2" type="123"/>
   <Atom name="HB3" type="123"/>
   <Atom name="CG" type="124"/>
   <Atom name="HG2" type="125"/>
   <Atom name="HG3" type="125"/>
   <Atom name="CD" type="126"/>
   <Atom name="OE1" type="127"/>
   <Atom name="OE2" type="128"/>
   <Atom name="C" type="129"/>
   <Atom name="O" type="130"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="13"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="7" to="9"/>
   <Bond from="7" to="10"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="12"/>
   <Bond from="13" to="14"/>
   <ExternalBond from="0"/>
   <ExternalBond from="13"/>
  </Residue>
  <Residue name="GLY">
   <Atom name="N" type="131"/>
   <Atom name="H" type="132"/>
   <Atom name="CA" type="133"/>
   <Atom name="HA2" type="134"/>
   <Atom name="HA3" type="134"/>
   <Atom name="C" type="135"/>
   <Atom name="O" type="136"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="5"/>
   <Bond from="5" to="6"/>
   <ExternalBond from="0"/>
   <ExternalBond from="5"/>
  </Residue>
  <Residue name="HID">
   <Atom name="N" type="137"/>
   <Atom name="H" type="138"/>
   <Atom name="CA" type="139"/>
   <Atom name="HA" type="140"/>
   <Atom name="CB" type="141"/>
   <Atom name="HB2" type="142"/>
   <Atom name="HB3" type="142"/>
   <Atom name="CG" type="143"/>
   <Atom name="ND1" type="144"/>
   <Atom name="HD1" type="145"/>
   <Atom name="CE1" type="146"/>
   <Atom name="HE1" type="147"/>
   <Atom name="NE2" type="148"/>
   <Atom name="CD2" type="149"/>
   <Atom name="HD2" type="150"/>
   <Atom name="C" type="151"/>
   <Atom name="O" type="152"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="15"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="7" to="13"/>
   <Bond from="8" to="9"/>
   <Bond from="8" to="10"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="12"/>
   <Bond from="12" to="13"/>
   <Bond from="13" to="14"/>
   <Bond from="15" to="16"/>
   <ExternalBond from="0"/>
   <ExternalBond from="15"/>
  </Residue>
  <Residue name="HIE">
   <Atom name="N" type="153"/>
   <Atom name="H" type="154"/>
   <Atom name="CA" type="155"/>
   <Atom name="HA" type="156"/>
   <Atom name="CB" type="157"/>
   <Atom name="HB2" type="158"/>
   <Atom name="HB3" type="158"/>
   <Atom name="CG" type="159"/>
   <Atom name="ND1" type="160"/>
   <Atom name="CE1" type="161"/>
   <Atom name="HE1" type="162"/>
   <Atom name="NE2" type="163"/>
   <Atom name="HE2" type="164"/>
   <Atom name="CD2" type="165"/>
   <Atom name="HD2" type="166"/>
   <Atom name="C" type="167"/>
   <Atom name="O" type="168"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="15"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="7" to="13"/>
   <Bond from="8" to="9"/>
   <Bond from="9" to="10"/>
   <Bond from="9" to="11"/>
   <Bond from="11" to="12"/>
   <Bond from="11" to="13"/>
   <Bond from="13" to="14"/>
   <Bond from="15" to="16"/>
   <ExternalBond from="0"/>
   <ExternalBond from="15"/>
  </Residue>
  <Residue name="HIP">
   <Atom name="N" type="169"/>
   <Atom name="H" type="170"/>
   <Atom name="CA" type="171"/>
   <Atom name="HA" type="172"/>
   <Atom name="CB" type="173"/>
   <Atom name="HB2" type="174"/>
   <Atom name="HB3" type="174"/>
   <Atom name="CG" type="175"/>
   <Atom name="ND1" type="176"/>
   <Atom name="HD1" type="177"/>
   <Atom name="CE1" type="178"/>
   <Atom name="HE1" type="179"/>
   <Atom name="NE2" type="180"/>
   <Atom name="HE2" type="181"/>
   <Atom name="CD2" type="182"/>
   <Atom name="HD2" type="183"/>
   <Atom name="C" type="184"/>
   <Atom name="O" type="185"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="16"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="7" to="14"/>
   <Bond from="8" to="9"/>
   <Bond from="8" to="10"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="12"/>
   <Bond from="12" to="13"/>
   <Bond from="12" to="14"/>
   <Bond from="14" to="15"/>
   <Bond from="16" to="17"/>
   <ExternalBond from="0"/>
   <ExternalBond from="16"/>
  </Residue>
  <Residue name="ILE">
   <Atom name="N" type="186"/>
   <Atom name="H" type="187"/>
   <Atom name="CA" type="188"/>
   <Atom name="HA" type="189"/>
   <Atom name="CB" type="190"/>
   <Atom name="HB" type="191"/>
   <Atom name="CG2" type="192"/>
   <Atom name="HG21" type="193"/>
   <Atom name="HG22" type="193"/>
   <Atom name="HG23" type="193"/>
   <Atom name="CG1" type="194"/>
   <Atom name="HG12" type="195"/>
   <Atom name="HG13" type="195"/>
   <Atom name="CD1" type="196"/>
   <Atom name="HD11" type="197"/>
   <Atom name="HD12" type="197"/>
   <Atom name="HD13" type="197"/>
   <Atom name="C" type="198"/>
   <Atom name="O" type="199"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="17"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="10"/>
   <Bond from="6" to="7"/>
   <Bond from="6" to="8"/>
   <Bond from="6" to="9"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="12"/>
   <Bond from="10" to="13"/>
   <Bond from="13" to="14"/>
   <Bond from="13" to="15"/>
   <Bond from="13" to="16"/>
   <Bond from="17" to="18"/>
   <ExternalBond from="0"/>
   <ExternalBond from="17"/>
  </Residue>
  <Residue name="K+">
   <Atom name="K+" type="1956"/>
  </Residue>
  <Residue name="LEU">
   <Atom name="N" type="200"/>
   <Atom name="H" type="201"/>
   <Atom name="CA" type="202"/>
   <Atom name="HA" type="203"/>
   <Atom name="CB" type="204"/>
   <Atom name="HB2" type="205"/>
   <Atom name="HB3" type="205"/>
   <Atom name="CG" type="206"/>
   <Atom name="HG" type="207"/>
   <Atom name="CD1" type="208"/>
   <Atom name="HD11" type="209"/>
   <Atom name="HD12" type="209"/>
   <Atom name="HD13" type="209"/>
   <Atom name="CD2" type="210"/>
   <Atom name="HD21" type="211"/>
   <Atom name="HD22" type="211"/>
   <Atom name="HD23" type="211"/>
   <Atom name="C" type="212"/>
   <Atom name="O" type="213"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="17"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="7" to="9"/>
   <Bond from="7" to="13"/>
   <Bond from="9" to="10"/>
   <Bond from="9" to="11"/>
   <Bond from="9" to="12"/>
   <Bond from="13" to="14"/>
   <Bond from="13" to="15"/>
   <Bond from="13" to="16"/>
   <Bond from="17" to="18"/>
   <ExternalBond from="0"/>
   <ExternalBond from="17"/>
  </Residue>
  <Residue name="LYN">
   <Atom name="N" type="214"/>
   <Atom name="H" type="215"/>
   <Atom name="CA" type="216"/>
   <Atom name="HA" type="217"/>
   <Atom name="CB" type="218"/>
   <Atom name="HB2" type="219"/>
   <Atom name="HB3" type="219"/>
   <Atom name="CG" type="220"/>
   <Atom name="HG2" type="221"/>
   <Atom name="HG3" type="221"/>
   <Atom name="CD" type="222"/>
   <Atom name="HD2" type="223"/>
   <Atom name="HD3" type="223"/>
   <Atom name="CE" type="224"/>
   <Atom name="HE2" type="225"/>
   <Atom name="HE3" type="225"/>
   <Atom name="NZ" type="226"/>
   <Atom name="HZ2" type="227"/>
   <Atom name="HZ3" type="227"/>
   <Atom name="C" type="228"/>
   <Atom name="O" type="229"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="19"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="7" to="9"/>
   <Bond from="7" to="10"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="12"/>
   <Bond from="10" to="13"/>
   <Bond from="13" to="14"/>
   <Bond from="13" to="15"/>
   <Bond from="13" to="16"/>
   <Bond from="16" to="17"/>
   <Bond from="16" to="18"/>
   <Bond from="19" to="20"/>
   <ExternalBond from="0"/>
   <ExternalBond from="19"/>
  </Residue>
  <Residue name="LYS">
   <Atom name="N" type="230"/>
   <Atom name="H" type="231"/>
   <Atom name="CA" type="232"/>
   <Atom name="HA" type="233"/>
   <Atom name="CB" type="234"/>
   <Atom name="HB2" type="235"/>
   <Atom name="HB3" type="235"/>
   <Atom name="CG" type="236"/>
   <Atom name="HG2" type="237"/>
   <Atom name="HG3" type="237"/>
   <Atom name="CD" type="238"/>
   <Atom name="HD2" type="239"/>
   <Atom name="HD3" type="239"/>
   <Atom name="CE" type="240"/>
   <Atom name="HE2" type="241"/>
   <Atom name="HE3" type="241"/>
   <Atom name="NZ" type="242"/>
   <Atom name="HZ1" type="243"/>
   <Atom name="HZ2" type="243"/>
   <Atom name="HZ3" type="243"/>
   <Atom name="C" type="244"/>
   <Atom name="O" type="245"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="20"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="7" to="9"/>
   <Bond from="7" to="10"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="12"/>
   <Bond from="10" to="13"/>
   <Bond from="13" to="14"/>
   <Bond from="13" to="15"/>
   <Bond from="13" to="16"/>
   <Bond from="16" to="17"/>
   <Bond from="16" to="18"/>
   <Bond from="16" to="19"/>
   <Bond from="20" to="21"/>
   <ExternalBond from="0"/>
   <ExternalBond from="20"/>
  </Residue>
  <Residue name="Li+">
   <Atom name="Li+" type="1957"/>
  </Residue>
  <Residue name="MET">
   <Atom name="N" type="246"/>
   <Atom name="H" type="247"/>
   <Atom name="CA" type="248"/>
   <Atom name="HA" type="249"/>
   <Atom name="CB" type="250"/>
   <Atom name="HB2" type="251"/>
   <Atom name="HB3" type="251"/>
   <Atom name="CG" type="252"/>
   <Atom name="HG2" type="253"/>
   <Atom name="HG3" type="253"/>
   <Atom name="SD" type="254"/>
   <Atom name="CE" type="255"/>
   <Atom name="HE1" type="256"/>
   <Atom name="HE2" type="256"/>
   <Atom name="HE3" type="256"/>
   <Atom name="C" type="257"/>
   <Atom name="O" type="258"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="15"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="7" to="9"/>
   <Bond from="7" to="10"/>
   <Bond from="10" to="11"/>
   <Bond from="11" to="12"/>
   <Bond from="11" to="13"/>
   <Bond from="11" to="14"/>
   <Bond from="15" to="16"/>
   <ExternalBond from="0"/>
   <ExternalBond from="15"/>
  </Residue>
  <Residue name="MG2">
   <Atom name="MG" type="1958"/>
  </Residue>
  <Residue name="NALA">
   <Atom name="N" type="714"/>
   <Atom name="H1" type="715"/>
   <Atom name="H2" type="715"/>
   <Atom name="H3" type="715"/>
   <Atom name="CA" type="716"/>
   <Atom name="HA" type="717"/>
   <Atom name="CB" type="718"/>
   <Atom name="HB1" type="719"/>
   <Atom name="HB2" type="719"/>
   <Atom name="HB3" type="719"/>
   <Atom name="C" type="720"/>
   <Atom name="O" type="721"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="0" to="3"/>
   <Bond from="0" to="4"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="10"/>
   <Bond from="6" to="7"/>
   <Bond from="6" to="8"/>
   <Bond from="6" to="9"/>
   <Bond from="10" to="11"/>
   <ExternalBond from="10"/>
  </Residue>
  <Residue name="NARG">
   <Atom name="N" type="722"/>
   <Atom name="H1" type="723"/>
   <Atom name="H2" type="723"/>
   <Atom name="H3" type="723"/>
   <Atom name="CA" type="724"/>
   <Atom name="HA" type="725"/>
   <Atom name="CB" type="726"/>
   <Atom name="HB2" type="727"/>
   <Atom name="HB3" type="727"/>
   <Atom name="CG" type="728"/>
   <Atom name="HG2" type="729"/>
   <Atom name="HG3" type="729"/>
   <Atom name="CD" type="730"/>
   <Atom name="HD2" type="731"/>
   <Atom name="HD3" type="731"/>
   <Atom name="NE" type="732"/>
   <Atom name="HE" type="733"/>
   <Atom name="CZ" type="734"/>
   <Atom name="NH1" type="735"/>
   <Atom name="HH11" type="736"/>
   <Atom name="HH12" type="736"/>
   <Atom name="NH2" type="737"/>
   <Atom name="HH21" type="738"/>
   <Atom name="HH22" type="738"/>
   <Atom name="C" type="739"/>
   <Atom name="O" type="740"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="0" to="3"/>
   <Bond from="0" to="4"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="24"/>
   <Bond from="6" to="7"/>
   <Bond from="6" to="8"/>
   <Bond from="6" to="9"/>
   <Bond from="9" to="10"/>
   <Bond from="9" to="11"/>
   <Bond from="9" to="12"/>
   <Bond from="12" to="13"/>
   <Bond from="12" to="14"/>
   <Bond from="12" to="15"/>
   <Bond from="15" to="16"/>
   <Bond from="15" to="17"/>
   <Bond from="17" to="18"/>
   <Bond from="17" to="21"/>
   <Bond from="18" to="19"/>
   <Bond from="18" to="20"/>
   <Bond from="21" to="22"/>
   <Bond from="21" to="23"/>
   <Bond from="24" to="25"/>
   <ExternalBond from="24"/>
  </Residue>
  <Residue name="NASN">
   <Atom name="N" type="741"/>
   <Atom name="H1" type="742"/>
   <Atom name="H2" type="742"/>
   <Atom name="H3" type="742"/>
   <Atom name="CA" type="743"/>
   <Atom name="HA" type="744"/>
   <Atom name="CB" type="745"/>
   <Atom name="HB2" type="746"/>
   <Atom name="HB3" type="746"/>
   <Atom name="CG" type="747"/>
   <Atom name="OD1" type="748"/>
   <Atom name="ND2" type="749"/>
   <Atom name="HD21" type="750"/>
   <Atom name="HD22" type="750"/>
   <Atom name="C" type="751"/>
   <Atom name="O" type="752"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="0" to="3"/>
   <Bond from="0" to="4"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="14"/>
   <Bond from="6" to="7"/>
   <Bond from="6" to="8"/>
   <Bond from="6" to="9"/>
   <Bond from="9" to="10"/>
   <Bond from="9" to="11"/>
   <Bond from="11" to="12"/>
   <Bond from="11" to="13"/>
   <Bond from="14" to="15"/>
   <ExternalBond from="14"/>
  </Residue>
  <Residue name="NASP">
   <Atom name="N" type="753"/>
   <Atom name="H1" type="754"/>
   <Atom name="H2" type="754"/>
   <Atom name="H3" type="754"/>
   <Atom name="CA" type="755"/>
   <Atom name="HA" type="756"/>
   <Atom name="CB" type="757"/>
   <Atom name="HB2" type="758"/>
   <Atom name="HB3" type="758"/>
   <Atom name="CG" type="759"/>
   <Atom name="OD1" type="760"/>
   <Atom name="OD2" type="761"/>
   <Atom name="C" type="762"/>
   <Atom name="O" type="763"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="0" to="3"/>
   <Bond from="0" to="4"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="12"/>
   <Bond from="6" to="7"/>
   <Bond from="6" to="8"/>
   <Bond from="6" to="9"/>
   <Bond from="9" to="10"/>
   <Bond from="9" to="11"/>
   <Bond from="12" to="13"/>
   <ExternalBond from="12"/>
  </Residue>
  <Residue name="NCYS">
   <Atom name="N" type="764"/>
   <Atom name="H1" type="765"/>
   <Atom name="H2" type="765"/>
   <Atom name="H3" type="765"/>
   <Atom name="CA" type="766"/>
   <Atom name="HA" type="767"/>
   <Atom name="CB" type="768"/>
   <Atom name="HB2" type="769"/>
   <Atom name="HB3" type="769"/>
   <Atom name="SG" type="770"/>
   <Atom name="HG" type="771"/>
   <Atom name="C" type="772"/>
   <Atom name="O" type="773"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="0" to="3"/>
   <Bond from="0" to="4"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="11"/>
   <Bond from="6" to="7"/>
   <Bond from="6" to="8"/>
   <Bond from="6" to="9"/>
   <Bond from="9" to="10"/>
   <Bond from="11" to="12"/>
   <ExternalBond from="11"/>
  </Residue>
  <Residue name="NCYX">
   <Atom name="N" type="774"/>
   <Atom name="H1" type="775"/>
   <Atom name="H2" type="775"/>
   <Atom name="H3" type="775"/>
   <Atom name="CA" type="776"/>
   <Atom name="HA" type="777"/>
   <Atom name="CB" type="778"/>
   <Atom name="HB2" type="779"/>
   <Atom name="HB3" type="779"/>
   <Atom name="SG" type="780"/>
   <Atom name="C" type="781"/>
   <Atom name="O" type="782"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="0" to="3"/>
   <Bond from="0" to="4"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="10"/>
   <Bond from="6" to="7"/>
   <Bond from="6" to="8"/>
   <Bond from="6" to="9"/>
   <Bond from="10" to="11"/>
   <ExternalBond from="10"/>
   <ExternalBond from="9"/>
  </Residue>
  <Residue name="NGLN">
   <Atom name="N" type="783"/>
   <Atom name="H1" type="784"/>
   <Atom name="H2" type="784"/>
   <Atom name="H3" type="784"/>
   <Atom name="CA" type="785"/>
   <Atom name="HA" type="786"/>
   <Atom name="CB" type="787"/>
   <Atom name="HB2" type="788"/>
   <Atom name="HB3" type="788"/>
   <Atom name="CG" type="789"/>
   <Atom name="HG2" type="790"/>
   <Atom name="HG3" type="790"/>
   <Atom name="CD" type="791"/>
   <Atom name="OE1" type="792"/>
   <Atom name="NE2" type="793"/>
   <Atom name="HE21" type="794"/>
   <Atom name="HE22" type="794"/>
   <Atom name="C" type="795"/>
   <Atom name="O" type="796"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="0" to="3"/>
   <Bond from="0" to="4"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="17"/>
   <Bond from="6" to="7"/>
   <Bond from="6" to="8"/>
   <Bond from="6" to="9"/>
   <Bond from="9" to="10"/>
   <Bond from="9" to="11"/>
   <Bond from="9" to="12"/>
   <Bond from="12" to="13"/>
   <Bond from="12" to="14"/>
   <Bond from="14" to="15"/>
   <Bond from="14" to="16"/>
   <Bond from="17" to="18"/>
   <ExternalBond from="17"/>
  </Residue>
  <Residue name="NGLU">
   <Atom name="N" type="797"/>
   <Atom name="H1" type="798"/>
   <Atom name="H2" type="798"/>
   <Atom name="H3" type="798"/>
   <Atom name="CA" type="799"/>
   <Atom name="HA" type="800"/>
   <Atom name="CB" type="801"/>
   <Atom name="HB2" type="802"/>
   <Atom name="HB3" type="802"/>
   <Atom name="CG" type="803"/>
   <Atom name="HG2" type="804"/>
   <Atom name="HG3" type="804"/>
   <Atom name="CD" type="805"/>
   <Atom name="OE1" type="806"/>
   <Atom name="OE2" type="807"/>
   <Atom name="C" type="808"/>
   <Atom name="O" type="809"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="0" to="3"/>
   <Bond from="0" to="4"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="15"/>
   <Bond from="6" to="7"/>
   <Bond from="6" to="8"/>
   <Bond from="6" to="9"/>
   <Bond from="9" to="10"/>
   <Bond from="9" to="11"/>
   <Bond from="9" to="12"/>
   <Bond from="12" to="13"/>
   <Bond from="12" to="14"/>
   <Bond from="15" to="16"/>
   <ExternalBond from="15"/>
  </Residue>
  <Residue name="NGLY">
   <Atom name="N" type="810"/>
   <Atom name="H1" type="811"/>
   <Atom name="H2" type="811"/>
   <Atom name="H3" type="811"/>
   <Atom name="CA" type="812"/>
   <Atom name="HA2" type="813"/>
   <Atom name="HA3" type="813"/>
   <Atom name="C" type="814"/>
   <Atom name="O" type="815"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="0" to="3"/>
   <Bond from="0" to="4"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <ExternalBond from="7"/>
  </Residue>
  <Residue name="NHE">
   <Atom name="N" type="704"/>
   <Atom name="HN1" type="705"/>
   <Atom name="HN2" type="705"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <ExternalBond from="0"/>
  </Residue>
  <Residue name="NHID">
   <Atom name="N" type="816"/>
   <Atom name="H1" type="817"/>
   <Atom name="H2" type="817"/>
   <Atom name="H3" type="817"/>
   <Atom name="CA" type="818"/>
   <Atom name="HA" type="819"/>
   <Atom name="CB" type="820"/>
   <Atom name="HB2" type="821"/>
   <Atom name="HB3" type="821"/>
   <Atom name="CG" type="822"/>
   <Atom name="ND1" type="823"/>
   <Atom name="HD1" type="824"/>
   <Atom name="CE1" type="825"/>
   <Atom name="HE1" type="826"/>
   <Atom name="NE2" type="827"/>
   <Atom name="CD2" type="828"/>
   <Atom name="HD2" type="829"/>
   <Atom name="C" type="830"/>
   <Atom name="O" type="831"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="0" to="3"/>
   <Bond from="0" to="4"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="17"/>
   <Bond from="6" to="7"/>
   <Bond from="6" to="8"/>
   <Bond from="6" to="9"/>
   <Bond from="9" to="10"/>
   <Bond from="9" to="15"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="12"/>
   <Bond from="12" to="13"/>
   <Bond from="12" to="14"/>
   <Bond from="14" to="15"/>
   <Bond from="15" to="16"/>
   <Bond from="17" to="18"/>
   <ExternalBond from="17"/>
  </Residue>
  <Residue name="NHIE">
   <Atom name="N" type="832"/>
   <Atom name="H1" type="833"/>
   <Atom name="H2" type="833"/>
   <Atom name="H3" type="833"/>
   <Atom name="CA" type="834"/>
   <Atom name="HA" type="835"/>
   <Atom name="CB" type="836"/>
   <Atom name="HB2" type="837"/>
   <Atom name="HB3" type="837"/>
   <Atom name="CG" type="838"/>
   <Atom name="ND1" type="839"/>
   <Atom name="CE1" type="840"/>
   <Atom name="HE1" type="841"/>
   <Atom name="NE2" type="842"/>
   <Atom name="HE2" type="843"/>
   <Atom name="CD2" type="844"/>
   <Atom name="HD2" type="845"/>
   <Atom name="C" type="846"/>
   <Atom name="O" type="847"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="0" to="3"/>
   <Bond from="0" to="4"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="17"/>
   <Bond from="6" to="7"/>
   <Bond from="6" to="8"/>
   <Bond from="6" to="9"/>
   <Bond from="9" to="10"/>
   <Bond from="9" to="15"/>
   <Bond from="10" to="11"/>
   <Bond from="11" to="12"/>
   <Bond from="11" to="13"/>
   <Bond from="13" to="14"/>
   <Bond from="13" to="15"/>
   <Bond from="15" to="16"/>
   <Bond from="17" to="18"/>
   <ExternalBond from="17"/>
  </Residue>
  <Residue name="NHIP">
   <Atom name="N" type="848"/>
   <Atom name="H1" type="849"/>
   <Atom name="H2" type="849"/>
   <Atom name="H3" type="849"/>
   <Atom name="CA" type="850"/>
   <Atom name="HA" type="851"/>
   <Atom name="CB" type="852"/>
   <Atom name="HB2" type="853"/>
   <Atom name="HB3" type="853"/>
   <Atom name="CG" type="854"/>
   <Atom name="ND1" type="855"/>
   <Atom name="HD1" type="856"/>
   <Atom name="CE1" type="857"/>
   <Atom name="HE1" type="858"/>
   <Atom name="NE2" type="859"/>
   <Atom name="HE2" type="860"/>
   <Atom name="CD2" type="861"/>
   <Atom name="HD2" type="862"/>
   <Atom name="C" type="863"/>
   <Atom name="O" type="864"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="0" to="3"/>
   <Bond from="0" to="4"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="18"/>
   <Bond from="6" to="7"/>
   <Bond from="6" to="8"/>
   <Bond from="6" to="9"/>
   <Bond from="9" to="10"/>
   <Bond from="9" to="16"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="12"/>
   <Bond from="12" to="13"/>
   <Bond from="12" to="14"/>
   <Bond from="14" to="15"/>
   <Bond from="14" to="16"/>
   <Bond from="16" to="17"/>
   <Bond from="18" to="19"/>
   <ExternalBond from="18"/>
  </Residue>
  <Residue name="NILE">
   <Atom name="N" type="865"/>
   <Atom name="H1" type="866"/>
   <Atom name="H2" type="866"/>
   <Atom name="H3" type="866"/>
   <Atom name="CA" type="867"/>
   <Atom name="HA" type="868"/>
   <Atom name="CB" type="869"/>
   <Atom name="HB" type="870"/>
   <Atom name="CG2" type="871"/>
   <Atom name="HG21" type="872"/>
   <Atom name="HG22" type="872"/>
   <Atom name="HG23" type="872"/>
   <Atom name="CG1" type="873"/>
   <Atom name="HG12" type="874"/>
   <Atom name="HG13" type="874"/>
   <Atom name="CD1" type="875"/>
   <Atom name="HD11" type="876"/>
   <Atom name="HD12" type="876"/>
   <Atom name="HD13" type="876"/>
   <Atom name="C" type="877"/>
   <Atom name="O" type="878"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="0" to="3"/>
   <Bond from="0" to="4"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="19"/>
   <Bond from="6" to="7"/>
   <Bond from="6" to="8"/>
   <Bond from="6" to="12"/>
   <Bond from="8" to="9"/>
   <Bond from="8" to="10"/>
   <Bond from="8" to="11"/>
   <Bond from="12" to="13"/>
   <Bond from="12" to="14"/>
   <Bond from="12" to="15"/>
   <Bond from="15" to="16"/>
   <Bond from="15" to="17"/>
   <Bond from="15" to="18"/>
   <Bond from="19" to="20"/>
   <ExternalBond from="19"/>
  </Residue>
  <Residue name="NLEU">
   <Atom name="N" type="879"/>
   <Atom name="H1" type="880"/>
   <Atom name="H2" type="880"/>
   <Atom name="H3" type="880"/>
   <Atom name="CA" type="881"/>
   <Atom name="HA" type="882"/>
   <Atom name="CB" type="883"/>
   <Atom name="HB2" type="884"/>
   <Atom name="HB3" type="884"/>
   <Atom name="CG" type="885"/>
   <Atom name="HG" type="886"/>
   <Atom name="CD1" type="887"/>
   <Atom name="HD11" type="888"/>
   <Atom name="HD12" type="888"/>
   <Atom name="HD13" type="888"/>
   <Atom name="CD2" type="889"/>
   <Atom name="HD21" type="890"/>
   <Atom name="HD22" type="890"/>
   <Atom name="HD23" type="890"/>
   <Atom name="C" type="891"/>
   <Atom name="O" type="892"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="0" to="3"/>
   <Bond from="0" to="4"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="19"/>
   <Bond from="6" to="7"/>
   <Bond from="6" to="8"/>
   <Bond from="6" to="9"/>
   <Bond from="9" to="10"/>
   <Bond from="9" to="11"/>
   <Bond from="9" to="15"/>
   <Bond from="11" to="12"/>
   <Bond from="11" to="13"/>
   <Bond from="11" to="14"/>
   <Bond from="15" to="16"/>
   <Bond from="15" to="17"/>
   <Bond from="15" to="18"/>
   <Bond from="19" to="20"/>
   <ExternalBond from="19"/>
  </Residue>
  <Residue name="NLYS">
   <Atom name="N" type="893"/>
   <Atom name="H1" type="894"/>
   <Atom name="H2" type="894"/>
   <Atom name="H3" type="894"/>
   <Atom name="CA" type="895"/>
   <Atom name="HA" type="896"/>
   <Atom name="CB" type="897"/>
   <Atom name="HB2" type="898"/>
   <Atom name="HB3" type="898"/>
   <Atom name="CG" type="899"/>
   <Atom name="HG2" type="900"/>
   <Atom name="HG3" type="900"/>
   <Atom name="CD" type="901"/>
   <Atom name="HD2" type="902"/>
   <Atom name="HD3" type="902"/>
   <Atom name="CE" type="903"/>
   <Atom name="HE2" type="904"/>
   <Atom name="HE3" type="904"/>
   <Atom name="NZ" type="905"/>
   <Atom name="HZ1" type="906"/>
   <Atom name="HZ2" type="906"/>
   <Atom name="HZ3" type="906"/>
   <Atom name="C" type="907"/>
   <Atom name="O" type="908"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="0" to="3"/>
   <Bond from="0" to="4"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="22"/>
   <Bond from="6" to="7"/>
   <Bond from="6" to="8"/>
   <Bond from="6" to="9"/>
   <Bond from="9" to="10"/>
   <Bond from="9" to="11"/>
   <Bond from="9" to="12"/>
   <Bond from="12" to="13"/>
   <Bond from="12" to="14"/>
   <Bond from="12" to="15"/>
   <Bond from="15" to="16"/>
   <Bond from="15" to="17"/>
   <Bond from="15" to="18"/>
   <Bond from="18" to="19"/>
   <Bond from="18" to="20"/>
   <Bond from="18" to="21"/>
   <Bond from="22" to="23"/>
   <ExternalBond from="22"/>
  </Residue>
  <Residue name="NME">
   <Atom name="N" type="706"/>
   <Atom name="H" type="707"/>
   <Atom name="CH3" type="708"/>
   <Atom name="HH31" type="709"/>
   <Atom name="HH32" type="709"/>
   <Atom name="HH33" type="709"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="5"/>
   <ExternalBond from="0"/>
  </Residue>
  <Residue name="NMET">
   <Atom name="N" type="909"/>
   <Atom name="H1" type="910"/>
   <Atom name="H2" type="910"/>
   <Atom name="H3" type="910"/>
   <Atom name="CA" type="911"/>
   <Atom name="HA" type="912"/>
   <Atom name="CB" type="913"/>
   <Atom name="HB2" type="914"/>
   <Atom name="HB3" type="914"/>
   <Atom name="CG" type="915"/>
   <Atom name="HG2" type="916"/>
   <Atom name="HG3" type="916"/>
   <Atom name="SD" type="917"/>
   <Atom name="CE" type="918"/>
   <Atom name="HE1" type="919"/>
   <Atom name="HE2" type="919"/>
   <Atom name="HE3" type="919"/>
   <Atom name="C" type="920"/>
   <Atom name="O" type="921"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="0" to="3"/>
   <Bond from="0" to="4"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="17"/>
   <Bond from="6" to="7"/>
   <Bond from="6" to="8"/>
   <Bond from="6" to="9"/>
   <Bond from="9" to="10"/>
   <Bond from="9" to="11"/>
   <Bond from="9" to="12"/>
   <Bond from="12" to="13"/>
   <Bond from="13" to="14"/>
   <Bond from="13" to="15"/>
   <Bond from="13" to="16"/>
   <Bond from="17" to="18"/>
   <ExternalBond from="17"/>
  </Residue>
  <Residue name="NPHE">
   <Atom name="N" type="922"/>
   <Atom name="H1" type="923"/>
   <Atom name="H2" type="923"/>
   <Atom name="H3" type="923"/>
   <Atom name="CA" type="924"/>
   <Atom name="HA" type="925"/>
   <Atom name="CB" type="926"/>
   <Atom name="HB2" type="927"/>
   <Atom name="HB3" type="927"/>
   <Atom name="CG" type="928"/>
   <Atom name="CD1" type="929"/>
   <Atom name="HD1" type="930"/>
   <Atom name="CE1" type="931"/>
   <Atom name="HE1" type="932"/>
   <Atom name="CZ" type="933"/>
   <Atom name="HZ" type="934"/>
   <Atom name="CE2" type="935"/>
   <Atom name="HE2" type="936"/>
   <Atom name="CD2" type="937"/>
   <Atom name="HD2" type="938"/>
   <Atom name="C" type="939"/>
   <Atom name="O" type="940"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="0" to="3"/>
   <Bond from="0" to="4"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="20"/>
   <Bond from="6" to="7"/>
   <Bond from="6" to="8"/>
   <Bond from="6" to="9"/>
   <Bond from="9" to="10"/>
   <Bond from="9" to="18"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="12"/>
   <Bond from="12" to="13"/>
   <Bond from="12" to="14"/>
   <Bond from="14" to="15"/>
   <Bond from="14" to="16"/>
   <Bond from="16" to="17"/>
   <Bond from="16" to="18"/>
   <Bond from="18" to="19"/>
   <Bond from="20" to="21"/>
   <ExternalBond from="20"/>
  </Residue>
  <Residue name="NPRO">
   <Atom name="N" type="941"/>
   <Atom name="H2" type="942"/>
   <Atom name="H3" type="942"/>
   <Atom name="CD" type="943"/>
   <Atom name="HD2" type="944"/>
   <Atom name="HD3" type="944"/>
   <Atom name="CG" type="945"/>
   <Atom name="HG2" type="946"/>
   <Atom name="HG3" type="946"/>
   <Atom name="CB" type="947"/>
   <Atom name="HB2" type="948"/>
   <Atom name="HB3" type="948"/>
   <Atom name="CA" type="949"/>
   <Atom name="HA" type="950"/>
   <Atom name="C" type="951"/>
   <Atom name="O" type="952"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="0" to="3"/>
   <Bond from="0" to="12"/>
   <Bond from="3" to="4"/>
   <Bond from="3" to="5"/>
   <Bond from="3" to="6"/>
   <Bond from="6" to="7"/>
   <Bond from="6" to="8"/>
   <Bond from="6" to="9"/>
   <Bond from="9" to="10"/>
   <Bond from="9" to="11"/>
   <Bond from="9" to="12"/>
   <Bond from="12" to="13"/>
   <Bond from="12" to="14"/>
   <Bond from="14" to="15"/>
   <ExternalBond from="14"/>
  </Residue>
  <Residue name="NSER">
   <Atom name="N" type="953"/>
   <Atom name="H1" type="954"/>
   <Atom name="H2" type="954"/>
   <Atom name="H3" type="954"/>
   <Atom name="CA" type="955"/>
   <Atom name="HA" type="956"/>
   <Atom name="CB" type="957"/>
   <Atom name="HB2" type="958"/>
   <Atom name="HB3" type="958"/>
   <Atom name="OG" type="959"/>
   <Atom name="HG" type="960"/>
   <Atom name="C" type="961"/>
   <Atom name="O" type="962"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="0" to="3"/>
   <Bond from="0" to="4"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="11"/>
   <Bond from="6" to="7"/>
   <Bond from="6" to="8"/>
   <Bond from="6" to="9"/>
   <Bond from="9" to="10"/>
   <Bond from="11" to="12"/>
   <ExternalBond from="11"/>
  </Residue>
  <Residue name="NTHR">
   <Atom name="N" type="963"/>
   <Atom name="H1" type="964"/>
   <Atom name="H2" type="964"/>
   <Atom name="H3" type="964"/>
   <Atom name="CA" type="965"/>
   <Atom name="HA" type="966"/>
   <Atom name="CB" type="967"/>
   <Atom name="HB" type="968"/>
   <Atom name="CG2" type="969"/>
   <Atom name="HG21" type="970"/>
   <Atom name="HG22" type="970"/>
   <Atom name="HG23" type="970"/>
   <Atom name="OG1" type="971"/>
   <Atom name="HG1" type="972"/>
   <Atom name="C" type="973"/>
   <Atom name="O" type="974"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="0" to="3"/>
   <Bond from="0" to="4"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="14"/>
   <Bond from="6" to="7"/>
   <Bond from="6" to="8"/>
   <Bond from="6" to="12"/>
   <Bond from="8" to="9"/>
   <Bond from="8" to="10"/>
   <Bond from="8" to="11"/>
   <Bond from="12" to="13"/>
   <Bond from="14" to="15"/>
   <ExternalBond from="14"/>
  </Residue>
  <Residue name="NTRP">
   <Atom name="N" type="975"/>
   <Atom name="H1" type="976"/>
   <Atom name="H2" type="976"/>
   <Atom name="H3" type="976"/>
   <Atom name="CA" type="977"/>
   <Atom name="HA" type="978"/>
   <Atom name="CB" type="979"/>
   <Atom name="HB2" type="980"/>
   <Atom name="HB3" type="980"/>
   <Atom name="CG" type="981"/>
   <Atom name="CD1" type="982"/>
   <Atom name="HD1" type="983"/>
   <Atom name="NE1" type="984"/>
   <Atom name="HE1" type="985"/>
   <Atom name="CE2" type="986"/>
   <Atom name="CZ2" type="987"/>
   <Atom name="HZ2" type="988"/>
   <Atom name="CH2" type="989"/>
   <Atom name="HH2" type="990"/>
   <Atom name="CZ3" type="991"/>
   <Atom name="HZ3" type="992"/>
   <Atom name="CE3" type="993"/>
   <Atom name="HE3" type="994"/>
   <Atom name="CD2" type="995"/>
   <Atom name="C" type="996"/>
   <Atom name="O" type="997"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="0" to="3"/>
   <Bond from="0" to="4"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="24"/>
   <Bond from="6" to="7"/>
   <Bond from="6" to="8"/>
   <Bond from="6" to="9"/>
   <Bond from="9" to="10"/>
   <Bond from="9" to="23"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="12"/>
   <Bond from="12" to="13"/>
   <Bond from="12" to="14"/>
   <Bond from="14" to="15"/>
   <Bond from="14" to="23"/>
   <Bond from="15" to="16"/>
   <Bond from="15" to="17"/>
   <Bond from="17" to="18"/>
   <Bond from="17" to="19"/>
   <Bond from="19" to="20"/>
   <Bond from="19" to="21"/>
   <Bond from="21" to="22"/>
   <Bond from="21" to="23"/>
   <Bond from="24" to="25"/>
   <ExternalBond from="24"/>
  </Residue>
  <Residue name="NTYR">
   <Atom name="N" type="998"/>
   <Atom name="H1" type="999"/>
   <Atom name="H2" type="999"/>
   <Atom name="H3" type="999"/>
   <Atom name="CA" type="1000"/>
   <Atom name="HA" type="1001"/>
   <Atom name="CB" type="1002"/>
   <Atom name="HB2" type="1003"/>
   <Atom name="HB3" type="1003"/>
   <Atom name="CG" type="1004"/>
   <Atom name="CD1" type="1005"/>
   <Atom name="HD1" type="1006"/>
   <Atom name="CE1" type="1007"/>
   <Atom name="HE1" type="1008"/>
   <Atom name="CZ" type="1009"/>
   <Atom name="OH" type="1010"/>
   <Atom name="HH" type="1011"/>
   <Atom name="CE2" type="1012"/>
   <Atom name="HE2" type="1013"/>
   <Atom name="CD2" type="1014"/>
   <Atom name="HD2" type="1015"/>
   <Atom name="C" type="1016"/>
   <Atom name="O" type="1017"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="0" to="3"/>
   <Bond from="0" to="4"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="21"/>
   <Bond from="6" to="7"/>
   <Bond from="6" to="8"/>
   <Bond from="6" to="9"/>
   <Bond from="9" to="10"/>
   <Bond from="9" to="19"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="12"/>
   <Bond from="12" to="13"/>
   <Bond from="12" to="14"/>
   <Bond from="14" to="15"/>
   <Bond from="14" to="17"/>
   <Bond from="15" to="16"/>
   <Bond from="17" to="18"/>
   <Bond from="17" to="19"/>
   <Bond from="19" to="20"/>
   <Bond from="21" to="22"/>
   <ExternalBond from="21"/>
  </Residue>
  <Residue name="NVAL">
   <Atom name="N" type="1018"/>
   <Atom name="H1" type="1019"/>
   <Atom name="H2" type="1019"/>
   <Atom name="H3" type="1019"/>
   <Atom name="CA" type="1020"/>
   <Atom name="HA" type="1021"/>
   <Atom name="CB" type="1022"/>
   <Atom name="HB" type="1023"/>
   <Atom name="CG1" type="1024"/>
   <Atom name="HG11" type="1025"/>
   <Atom name="HG12" type="1025"/>
   <Atom name="HG13" type="1025"/>
   <Atom name="CG2" type="1026"/>
   <Atom name="HG21" type="1027"/>
   <Atom name="HG22" type="1027"/>
   <Atom name="HG23" type="1027"/>
   <Atom name="C" type="1028"/>
   <Atom name="O" type="1029"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="0" to="3"/>
   <Bond from="0" to="4"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="16"/>
   <Bond from="6" to="7"/>
   <Bond from="6" to="8"/>
   <Bond from="6" to="12"/>
   <Bond from="8" to="9"/>
   <Bond from="8" to="10"/>
   <Bond from="8" to="11"/>
   <Bond from="12" to="13"/>
   <Bond from="12" to="14"/>
   <Bond from="12" to="15"/>
   <Bond from="16" to="17"/>
   <ExternalBond from="16"/>
  </Residue>
  <Residue name="Na+">
   <Atom name="Na+" type="1959"/>
  </Residue>
  <Residue name="PHE">
   <Atom name="N" type="259"/>
   <Atom name="H" type="260"/>
   <Atom name="CA" type="261"/>
   <Atom name="HA" type="262"/>
   <Atom name="CB" type="263"/>
   <Atom name="HB2" type="264"/>
   <Atom name="HB3" type="264"/>
   <Atom name="CG" type="265"/>
   <Atom name="CD1" type="266"/>
   <Atom name="HD1" type="267"/>
   <Atom name="CE1" type="268"/>
   <Atom name="HE1" type="269"/>
   <Atom name="CZ" type="270"/>
   <Atom name="HZ" type="271"/>
   <Atom name="CE2" type="272"/>
   <Atom name="HE2" type="273"/>
   <Atom name="CD2" type="274"/>
   <Atom name="HD2" type="275"/>
   <Atom name="C" type="276"/>
   <Atom name="O" type="277"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="18"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="7" to="16"/>
   <Bond from="8" to="9"/>
   <Bond from="8" to="10"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="12"/>
   <Bond from="12" to="13"/>
   <Bond from="12" to="14"/>
   <Bond from="14" to="15"/>
   <Bond from="14" to="16"/>
   <Bond from="16" to="17"/>
   <Bond from="18" to="19"/>
   <ExternalBond from="0"/>
   <ExternalBond from="18"/>
  </Residue>
  <Residue name="PRO">
   <Atom name="N" type="278"/>
   <Atom name="CD" type="279"/>
   <Atom name="HD2" type="280"/>
   <Atom name="HD3" type="280"/>
   <Atom name="CG" type="281"/>
   <Atom name="HG2" type="282"/>
   <Atom name="HG3" type="282"/>
   <Atom name="CB" type="283"/>
   <Atom name="HB2" type="284"/>
   <Atom name="HB3" type="284"/>
   <Atom name="CA" type="285"/>
   <Atom name="HA" type="286"/>
   <Atom name="C" type="287"/>
   <Atom name="O" type="288"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="10"/>
   <Bond from="1" to="2"/>
   <Bond from="1" to="3"/>
   <Bond from="1" to="4"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="7" to="9"/>
   <Bond from="7" to="10"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="12"/>
   <Bond from="12" to="13"/>
   <ExternalBond from="0"/>
   <ExternalBond from="12"/>
  </Residue>
  <Residue name="RA">
   <Atom name="P" type="1478"/>
   <Atom name="O1P" type="1479"/>
   <Atom name="O2P" type="1480"/>
   <Atom name="O5'" type="1481"/>
   <Atom name="C5'" type="1482"/>
   <Atom name="H5'1" type="1483"/>
   <Atom name="H5'2" type="1483"/>
   <Atom name="C4'" type="1484"/>
   <Atom name="H4'" type="1485"/>
   <Atom name="O4'" type="1486"/>
   <Atom name="C1'" type="1487"/>
   <Atom name="H1'" type="1488"/>
   <Atom name="N9" type="1489"/>
   <Atom name="C8" type="1490"/>
   <Atom name="H8" type="1491"/>
   <Atom name="N7" type="1492"/>
   <Atom name="C5" type="1493"/>
   <Atom name="C6" type="1494"/>
   <Atom name="N6" type="1495"/>
   <Atom name="H61" type="1496"/>
   <Atom name="H62" type="1496"/>
   <Atom name="N1" type="1497"/>
   <Atom name="C2" type="1498"/>
   <Atom name="H2" type="1499"/>
   <Atom name="N3" type="1500"/>
   <Atom name="C4" type="1501"/>
   <Atom name="C3'" type="1502"/>
   <Atom name="H3'" type="1503"/>
   <Atom name="C2'" type="1504"/>
   <Atom name="H2'1" type="1505"/>
   <Atom name="O2'" type="1506"/>
   <Atom name="HO'2" type="1507"/>
   <Atom name="O3'" type="1508"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="0" to="3"/>
   <Bond from="3" to="4"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="7" to="9"/>
   <Bond from="7" to="26"/>
   <Bond from="9" to="10"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="12"/>
   <Bond from="10" to="28"/>
   <Bond from="12" to="13"/>
   <Bond from="12" to="25"/>
   <Bond from="13" to="14"/>
   <Bond from="13" to="15"/>
   <Bond from="15" to="16"/>
   <Bond from="16" to="17"/>
   <Bond from="16" to="25"/>
   <Bond from="17" to="18"/>
   <Bond from="17" to="21"/>
   <Bond from="18" to="19"/>
   <Bond from="18" to="20"/>
   <Bond from="21" to="22"/>
   <Bond from="22" to="23"/>
   <Bond from="22" to="24"/>
   <Bond from="24" to="25"/>
   <Bond from="26" to="27"/>
   <Bond from="26" to="28"/>
   <Bond from="26" to="32"/>
   <Bond from="28" to="29"/>
   <Bond from="28" to="30"/>
   <Bond from="30" to="31"/>
   <ExternalBond from="0"/>
   <ExternalBond from="32"/>
  </Residue>
  <Residue name="RA3">
   <Atom name="P" type="1509"/>
   <Atom name="O1P" type="1510"/>
   <Atom name="O2P" type="1511"/>
   <Atom name="O5'" type="1512"/>
   <Atom name="C5'" type="1513"/>
   <Atom name="H5'1" type="1514"/>
   <Atom name="H5'2" type="1514"/>
   <Atom name="C4'" type="1515"/>
   <Atom name="H4'" type="1516"/>
   <Atom name="O4'" type="1517"/>
   <Atom name="C1'" type="1518"/>
   <Atom name="H1'" type="1519"/>
   <Atom name="N9" type="1520"/>
   <Atom name="C8" type="1521"/>
   <Atom name="H8" type="1522"/>
   <Atom name="N7" type="1523"/>
   <Atom name="C5" type="1524"/>
   <Atom name="C6" type="1525"/>
   <Atom name="N6" type="1526"/>
   <Atom name="H61" type="1527"/>
   <Atom name="H62" type="1527"/>
   <Atom name="N1" type="1528"/>
   <Atom name="C2" type="1529"/>
   <Atom name="H2" type="1530"/>
   <Atom name="N3" type="1531"/>
   <Atom name="C4" type="1532"/>
   <Atom name="C3'" type="1533"/>
   <Atom name="H3'" type="1534"/>
   <Atom name="C2'" type="1535"/>
   <Atom name="H2'1" type="1536"/>
   <Atom name="O2'" type="1537"/>
   <Atom name="HO'2" type="1538"/>
   <Atom name="O3'" type="1539"/>
   <Atom name="H3T" type="1540"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="0" to="3"/>
   <Bond from="3" to="4"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="7" to="9"/>
   <Bond from="7" to="26"/>
   <Bond from="9" to="10"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="12"/>
   <Bond from="10" to="28"/>
   <Bond from="12" to="13"/>
   <Bond from="12" to="25"/>
   <Bond from="13" to="14"/>
   <Bond from="13" to="15"/>
   <Bond from="15" to="16"/>
   <Bond from="16" to="17"/>
   <Bond from="16" to="25"/>
   <Bond from="17" to="18"/>
   <Bond from="17" to="21"/>
   <Bond from="18" to="19"/>
   <Bond from="18" to="20"/>
   <Bond from="21" to="22"/>
   <Bond from="22" to="23"/>
   <Bond from="22" to="24"/>
   <Bond from="24" to="25"/>
   <Bond from="26" to="27"/>
   <Bond from="26" to="28"/>
   <Bond from="26" to="32"/>
   <Bond from="28" to="29"/>
   <Bond from="28" to="30"/>
   <Bond from="30" to="31"/>
   <Bond from="32" to="33"/>
   <ExternalBond from="0"/>
  </Residue>
  <Residue name="RA5">
   <Atom name="H5T" type="1541"/>
   <Atom name="O5'" type="1542"/>
   <Atom name="C5'" type="1543"/>
   <Atom name="H5'1" type="1544"/>
   <Atom name="H5'2" type="1544"/>
   <Atom name="C4'" type="1545"/>
   <Atom name="H4'" type="1546"/>
   <Atom name="O4'" type="1547"/>
   <Atom name="C1'" type="1548"/>
   <Atom name="H1'" type="1549"/>
   <Atom name="N9" type="1550"/>
   <Atom name="C8" type="1551"/>
   <Atom name="H8" type="1552"/>
   <Atom name="N7" type="1553"/>
   <Atom name="C5" type="1554"/>
   <Atom name="C6" type="1555"/>
   <Atom name="N6" type="1556"/>
   <Atom name="H61" type="1557"/>
   <Atom name="H62" type="1557"/>
   <Atom name="N1" type="1558"/>
   <Atom name="C2" type="1559"/>
   <Atom name="H2" type="1560"/>
   <Atom name="N3" type="1561"/>
   <Atom name="C4" type="1562"/>
   <Atom name="C3'" type="1563"/>
   <Atom name="H3'" type="1564"/>
   <Atom name="C2'" type="1565"/>
   <Atom name="H2'1" type="1566"/>
   <Atom name="O2'" type="1567"/>
   <Atom name="HO'2" type="1568"/>
   <Atom name="O3'" type="1569"/>
   <Bond from="0" to="1"/>
   <Bond from="1" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="5"/>
   <Bond from="5" to="6"/>
   <Bond from="5" to="7"/>
   <Bond from="5" to="24"/>
   <Bond from="7" to="8"/>
   <Bond from="8" to="9"/>
   <Bond from="8" to="10"/>
   <Bond from="8" to="26"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="23"/>
   <Bond from="11" to="12"/>
   <Bond from="11" to="13"/>
   <Bond from="13" to="14"/>
   <Bond from="14" to="15"/>
   <Bond from="14" to="23"/>
   <Bond from="15" to="16"/>
   <Bond from="15" to="19"/>
   <Bond from="16" to="17"/>
   <Bond from="16" to="18"/>
   <Bond from="19" to="20"/>
   <Bond from="20" to="21"/>
   <Bond from="20" to="22"/>
   <Bond from="22" to="23"/>
   <Bond from="24" to="25"/>
   <Bond from="24" to="26"/>
   <Bond from="24" to="30"/>
   <Bond from="26" to="27"/>
   <Bond from="26" to="28"/>
   <Bond from="28" to="29"/>
   <ExternalBond from="30"/>
  </Residue>
  <Residue name="RAN">
   <Atom name="H5T" type="1570"/>
   <Atom name="O5'" type="1571"/>
   <Atom name="C5'" type="1572"/>
   <Atom name="H5'1" type="1573"/>
   <Atom name="H5'2" type="1573"/>
   <Atom name="C4'" type="1574"/>
   <Atom name="H4'" type="1575"/>
   <Atom name="O4'" type="1576"/>
   <Atom name="C1'" type="1577"/>
   <Atom name="H1'" type="1578"/>
   <Atom name="N9" type="1579"/>
   <Atom name="C8" type="1580"/>
   <Atom name="H8" type="1581"/>
   <Atom name="N7" type="1582"/>
   <Atom name="C5" type="1583"/>
   <Atom name="C6" type="1584"/>
   <Atom name="N6" type="1585"/>
   <Atom name="H61" type="1586"/>
   <Atom name="H62" type="1586"/>
   <Atom name="N1" type="1587"/>
   <Atom name="C2" type="1588"/>
   <Atom name="H2" type="1589"/>
   <Atom name="N3" type="1590"/>
   <Atom name="C4" type="1591"/>
   <Atom name="C3'" type="1592"/>
   <Atom name="H3'" type="1593"/>
   <Atom name="C2'" type="1594"/>
   <Atom name="H2'1" type="1595"/>
   <Atom name="O2'" type="1596"/>
   <Atom name="HO'2" type="1597"/>
   <Atom name="O3'" type="1598"/>
   <Atom name="H3T" type="1599"/>
   <Bond from="0" to="1"/>
   <Bond from="1" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="5"/>
   <Bond from="5" to="6"/>
   <Bond from="5" to="7"/>
   <Bond from="5" to="24"/>
   <Bond from="7" to="8"/>
   <Bond from="8" to="9"/>
   <Bond from="8" to="10"/>
   <Bond from="8" to="26"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="23"/>
   <Bond from="11" to="12"/>
   <Bond from="11" to="13"/>
   <Bond from="13" to="14"/>
   <Bond from="14" to="15"/>
   <Bond from="14" to="23"/>
   <Bond from="15" to="16"/>
   <Bond from="15" to="19"/>
   <Bond from="16" to="17"/>
   <Bond from="16" to="18"/>
   <Bond from="19" to="20"/>
   <Bond from="20" to="21"/>
   <Bond from="20" to="22"/>
   <Bond from="22" to="23"/>
   <Bond from="24" to="25"/>
   <Bond from="24" to="26"/>
   <Bond from="24" to="30"/>
   <Bond from="26" to="27"/>
   <Bond from="26" to="28"/>
   <Bond from="28" to="29"/>
   <Bond from="30" to="31"/>
  </Residue>
  <Residue name="RC">
   <Atom name="P" type="1600"/>
   <Atom name="O1P" type="1601"/>
   <Atom name="O2P" type="1602"/>
   <Atom name="O5'" type="1603"/>
   <Atom name="C5'" type="1604"/>
   <Atom name="H5'1" type="1605"/>
   <Atom name="H5'2" type="1605"/>
   <Atom name="C4'" type="1606"/>
   <Atom name="H4'" type="1607"/>
   <Atom name="O4'" type="1608"/>
   <Atom name="C1'" type="1609"/>
   <Atom name="H1'" type="1610"/>
   <Atom name="N1" type="1611"/>
   <Atom name="C6" type="1612"/>
   <Atom name="H6" type="1613"/>
   <Atom name="C5" type="1614"/>
   <Atom name="H5" type="1615"/>
   <Atom name="C4" type="1616"/>
   <Atom name="N4" type="1617"/>
   <Atom name="H41" type="1618"/>
   <Atom name="H42" type="1618"/>
   <Atom name="N3" type="1619"/>
   <Atom name="C2" type="1620"/>
   <Atom name="O2" type="1621"/>
   <Atom name="C3'" type="1622"/>
   <Atom name="H3'" type="1623"/>
   <Atom name="C2'" type="1624"/>
   <Atom name="H2'1" type="1625"/>
   <Atom name="O2'" type="1626"/>
   <Atom name="HO'2" type="1627"/>
   <Atom name="O3'" type="1628"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="0" to="3"/>
   <Bond from="3" to="4"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="7" to="9"/>
   <Bond from="7" to="24"/>
   <Bond from="9" to="10"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="12"/>
   <Bond from="10" to="26"/>
   <Bond from="12" to="13"/>
   <Bond from="12" to="22"/>
   <Bond from="13" to="14"/>
   <Bond from="13" to="15"/>
   <Bond from="15" to="16"/>
   <Bond from="15" to="17"/>
   <Bond from="17" to="18"/>
   <Bond from="17" to="21"/>
   <Bond from="18" to="19"/>
   <Bond from="18" to="20"/>
   <Bond from="21" to="22"/>
   <Bond from="22" to="23"/>
   <Bond from="24" to="25"/>
   <Bond from="24" to="26"/>
   <Bond from="24" to="30"/>
   <Bond from="26" to="27"/>
   <Bond from="26" to="28"/>
   <Bond from="28" to="29"/>
   <ExternalBond from="0"/>
   <ExternalBond from="30"/>
  </Residue>
  <Residue name="RC3">
   <Atom name="P" type="1629"/>
   <Atom name="O1P" type="1630"/>
   <Atom name="O2P" type="1631"/>
   <Atom name="O5'" type="1632"/>
   <Atom name="C5'" type="1633"/>
   <Atom name="H5'1" type="1634"/>
   <Atom name="H5'2" type="1634"/>
   <Atom name="C4'" type="1635"/>
   <Atom name="H4'" type="1636"/>
   <Atom name="O4'" type="1637"/>
   <Atom name="C1'" type="1638"/>
   <Atom name="H1'" type="1639"/>
   <Atom name="N1" type="1640"/>
   <Atom name="C6" type="1641"/>
   <Atom name="H6" type="1642"/>
   <Atom name="C5" type="1643"/>
   <Atom name="H5" type="1644"/>
   <Atom name="C4" type="1645"/>
   <Atom name="N4" type="1646"/>
   <Atom name="H41" type="1647"/>
   <Atom name="H42" type="1647"/>
   <Atom name="N3" type="1648"/>
   <Atom name="C2" type="1649"/>
   <Atom name="O2" type="1650"/>
   <Atom name="C3'" type="1651"/>
   <Atom name="H3'" type="1652"/>
   <Atom name="C2'" type="1653"/>
   <Atom name="H2'1" type="1654"/>
   <Atom name="O2'" type="1655"/>
   <Atom name="HO'2" type="1656"/>
   <Atom name="O3'" type="1657"/>
   <Atom name="H3T" type="1658"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="0" to="3"/>
   <Bond from="3" to="4"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="7" to="9"/>
   <Bond from="7" to="24"/>
   <Bond from="9" to="10"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="12"/>
   <Bond from="10" to="26"/>
   <Bond from="12" to="13"/>
   <Bond from="12" to="22"/>
   <Bond from="13" to="14"/>
   <Bond from="13" to="15"/>
   <Bond from="15" to="16"/>
   <Bond from="15" to="17"/>
   <Bond from="17" to="18"/>
   <Bond from="17" to="21"/>
   <Bond from="18" to="19"/>
   <Bond from="18" to="20"/>
   <Bond from="21" to="22"/>
   <Bond from="22" to="23"/>
   <Bond from="24" to="25"/>
   <Bond from="24" to="26"/>
   <Bond from="24" to="30"/>
   <Bond from="26" to="27"/>
   <Bond from="26" to="28"/>
   <Bond from="28" to="29"/>
   <Bond from="30" to="31"/>
   <ExternalBond from="0"/>
  </Residue>
  <Residue name="RC5">
   <Atom name="H5T" type="1659"/>
   <Atom name="O5'" type="1660"/>
   <Atom name="C5'" type="1661"/>
   <Atom name="H5'1" type="1662"/>
   <Atom name="H5'2" type="1662"/>
   <Atom name="C4'" type="1663"/>
   <Atom name="H4'" type="1664"/>
   <Atom name="O4'" type="1665"/>
   <Atom name="C1'" type="1666"/>
   <Atom name="H1'" type="1667"/>
   <Atom name="N1" type="1668"/>
   <Atom name="C6" type="1669"/>
   <Atom name="H6" type="1670"/>
   <Atom name="C5" type="1671"/>
   <Atom name="H5" type="1672"/>
   <Atom name="C4" type="1673"/>
   <Atom name="N4" type="1674"/>
   <Atom name="H41" type="1675"/>
   <Atom name="H42" type="1675"/>
   <Atom name="N3" type="1676"/>
   <Atom name="C2" type="1677"/>
   <Atom name="O2" type="1678"/>
   <Atom name="C3'" type="1679"/>
   <Atom name="H3'" type="1680"/>
   <Atom name="C2'" type="1681"/>
   <Atom name="H2'1" type="1682"/>
   <Atom name="O2'" type="1683"/>
   <Atom name="HO'2" type="1684"/>
   <Atom name="O3'" type="1685"/>
   <Bond from="0" to="1"/>
   <Bond from="1" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="5"/>
   <Bond from="5" to="6"/>
   <Bond from="5" to="7"/>
   <Bond from="5" to="22"/>
   <Bond from="7" to="8"/>
   <Bond from="8" to="9"/>
   <Bond from="8" to="10"/>
   <Bond from="8" to="24"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="20"/>
   <Bond from="11" to="12"/>
   <Bond from="11" to="13"/>
   <Bond from="13" to="14"/>
   <Bond from="13" to="15"/>
   <Bond from="15" to="16"/>
   <Bond from="15" to="19"/>
   <Bond from="16" to="17"/>
   <Bond from="16" to="18"/>
   <Bond from="19" to="20"/>
   <Bond from="20" to="21"/>
   <Bond from="22" to="23"/>
   <Bond from="22" to="24"/>
   <Bond from="22" to="28"/>
   <Bond from="24" to="25"/>
   <Bond from="24" to="26"/>
   <Bond from="26" to="27"/>
   <ExternalBond from="28"/>
  </Residue>
  <Residue name="RCN">
   <Atom name="H5T" type="1686"/>
   <Atom name="O5'" type="1687"/>
   <Atom name="C5'" type="1688"/>
   <Atom name="H5'1" type="1689"/>
   <Atom name="H5'2" type="1689"/>
   <Atom name="C4'" type="1690"/>
   <Atom name="H4'" type="1691"/>
   <Atom name="O4'" type="1692"/>
   <Atom name="C1'" type="1693"/>
   <Atom name="H1'" type="1694"/>
   <Atom name="N1" type="1695"/>
   <Atom name="C6" type="1696"/>
   <Atom name="H6" type="1697"/>
   <Atom name="C5" type="1698"/>
   <Atom name="H5" type="1699"/>
   <Atom name="C4" type="1700"/>
   <Atom name="N4" type="1701"/>
   <Atom name="H41" type="1702"/>
   <Atom name="H42" type="1702"/>
   <Atom name="N3" type="1703"/>
   <Atom name="C2" type="1704"/>
   <Atom name="O2" type="1705"/>
   <Atom name="C3'" type="1706"/>
   <Atom name="H3'" type="1707"/>
   <Atom name="C2'" type="1708"/>
   <Atom name="H2'1" type="1709"/>
   <Atom name="O2'" type="1710"/>
   <Atom name="HO'2" type="1711"/>
   <Atom name="O3'" type="1712"/>
   <Atom name="H3T" type="1713"/>
   <Bond from="0" to="1"/>
   <Bond from="1" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="5"/>
   <Bond from="5" to="6"/>
   <Bond from="5" to="7"/>
   <Bond from="5" to="22"/>
   <Bond from="7" to="8"/>
   <Bond from="8" to="9"/>
   <Bond from="8" to="10"/>
   <Bond from="8" to="24"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="20"/>
   <Bond from="11" to="12"/>
   <Bond from="11" to="13"/>
   <Bond from="13" to="14"/>
   <Bond from="13" to="15"/>
   <Bond from="15" to="16"/>
   <Bond from="15" to="19"/>
   <Bond from="16" to="17"/>
   <Bond from="16" to="18"/>
   <Bond from="19" to="20"/>
   <Bond from="20" to="21"/>
   <Bond from="22" to="23"/>
   <Bond from="22" to="24"/>
   <Bond from="22" to="28"/>
   <Bond from="24" to="25"/>
   <Bond from="24" to="26"/>
   <Bond from="26" to="27"/>
   <Bond from="28" to="29"/>
  </Residue>
  <Residue name="RG">
   <Atom name="P" type="1714"/>
   <Atom name="O1P" type="1715"/>
   <Atom name="O2P" type="1716"/>
   <Atom name="O5'" type="1717"/>
   <Atom name="C5'" type="1718"/>
   <Atom name="H5'1" type="1719"/>
   <Atom name="H5'2" type="1719"/>
   <Atom name="C4'" type="1720"/>
   <Atom name="H4'" type="1721"/>
   <Atom name="O4'" type="1722"/>
   <Atom name="C1'" type="1723"/>
   <Atom name="H1'" type="1724"/>
   <Atom name="N9" type="1725"/>
   <Atom name="C8" type="1726"/>
   <Atom name="H8" type="1727"/>
   <Atom name="N7" type="1728"/>
   <Atom name="C5" type="1729"/>
   <Atom name="C6" type="1730"/>
   <Atom name="O6" type="1731"/>
   <Atom name="N1" type="1732"/>
   <Atom name="H1" type="1733"/>
   <Atom name="C2" type="1734"/>
   <Atom name="N2" type="1735"/>
   <Atom name="H21" type="1736"/>
   <Atom name="H22" type="1736"/>
   <Atom name="N3" type="1737"/>
   <Atom name="C4" type="1738"/>
   <Atom name="C3'" type="1739"/>
   <Atom name="H3'" type="1740"/>
   <Atom name="C2'" type="1741"/>
   <Atom name="H2'1" type="1742"/>
   <Atom name="O2'" type="1743"/>
   <Atom name="HO'2" type="1744"/>
   <Atom name="O3'" type="1745"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="0" to="3"/>
   <Bond from="3" to="4"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="7" to="9"/>
   <Bond from="7" to="27"/>
   <Bond from="9" to="10"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="12"/>
   <Bond from="10" to="29"/>
   <Bond from="12" to="13"/>
   <Bond from="12" to="26"/>
   <Bond from="13" to="14"/>
   <Bond from="13" to="15"/>
   <Bond from="15" to="16"/>
   <Bond from="16" to="17"/>
   <Bond from="16" to="26"/>
   <Bond from="17" to="18"/>
   <Bond from="17" to="19"/>
   <Bond from="19" to="20"/>
   <Bond from="19" to="21"/>
   <Bond from="21" to="22"/>
   <Bond from="21" to="25"/>
   <Bond from="22" to="23"/>
   <Bond from="22" to="24"/>
   <Bond from="25" to="26"/>
   <Bond from="27" to="28"/>
   <Bond from="27" to="29"/>
   <Bond from="27" to="33"/>
   <Bond from="29" to="30"/>
   <Bond from="29" to="31"/>
   <Bond from="31" to="32"/>
   <ExternalBond from="0"/>
   <ExternalBond from="33"/>
  </Residue>
  <Residue name="RG3">
   <Atom name="P" type="1746"/>
   <Atom name="O1P" type="1747"/>
   <Atom name="O2P" type="1748"/>
   <Atom name="O5'" type="1749"/>
   <Atom name="C5'" type="1750"/>
   <Atom name="H5'1" type="1751"/>
   <Atom name="H5'2" type="1751"/>
   <Atom name="C4'" type="1752"/>
   <Atom name="H4'" type="1753"/>
   <Atom name="O4'" type="1754"/>
   <Atom name="C1'" type="1755"/>
   <Atom name="H1'" type="1756"/>
   <Atom name="N9" type="1757"/>
   <Atom name="C8" type="1758"/>
   <Atom name="H8" type="1759"/>
   <Atom name="N7" type="1760"/>
   <Atom name="C5" type="1761"/>
   <Atom name="C6" type="1762"/>
   <Atom name="O6" type="1763"/>
   <Atom name="N1" type="1764"/>
   <Atom name="H1" type="1765"/>
   <Atom name="C2" type="1766"/>
   <Atom name="N2" type="1767"/>
   <Atom name="H21" type="1768"/>
   <Atom name="H22" type="1768"/>
   <Atom name="N3" type="1769"/>
   <Atom name="C4" type="1770"/>
   <Atom name="C3'" type="1771"/>
   <Atom name="H3'" type="1772"/>
   <Atom name="C2'" type="1773"/>
   <Atom name="H2'1" type="1774"/>
   <Atom name="O2'" type="1775"/>
   <Atom name="HO'2" type="1776"/>
   <Atom name="O3'" type="1777"/>
   <Atom name="H3T" type="1778"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="0" to="3"/>
   <Bond from="3" to="4"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="7" to="9"/>
   <Bond from="7" to="27"/>
   <Bond from="9" to="10"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="12"/>
   <Bond from="10" to="29"/>
   <Bond from="12" to="13"/>
   <Bond from="12" to="26"/>
   <Bond from="13" to="14"/>
   <Bond from="13" to="15"/>
   <Bond from="15" to="16"/>
   <Bond from="16" to="17"/>
   <Bond from="16" to="26"/>
   <Bond from="17" to="18"/>
   <Bond from="17" to="19"/>
   <Bond from="19" to="20"/>
   <Bond from="19" to="21"/>
   <Bond from="21" to="22"/>
   <Bond from="21" to="25"/>
   <Bond from="22" to="23"/>
   <Bond from="22" to="24"/>
   <Bond from="25" to="26"/>
   <Bond from="27" to="28"/>
   <Bond from="27" to="29"/>
   <Bond from="27" to="33"/>
   <Bond from="29" to="30"/>
   <Bond from="29" to="31"/>
   <Bond from="31" to="32"/>
   <Bond from="33" to="34"/>
   <ExternalBond from="0"/>
  </Residue>
  <Residue name="RG5">
   <Atom name="H5T" type="1779"/>
   <Atom name="O5'" type="1780"/>
   <Atom name="C5'" type="1781"/>
   <Atom name="H5'1" type="1782"/>
   <Atom name="H5'2" type="1782"/>
   <Atom name="C4'" type="1783"/>
   <Atom name="H4'" type="1784"/>
   <Atom name="O4'" type="1785"/>
   <Atom name="C1'" type="1786"/>
   <Atom name="H1'" type="1787"/>
   <Atom name="N9" type="1788"/>
   <Atom name="C8" type="1789"/>
   <Atom name="H8" type="1790"/>
   <Atom name="N7" type="1791"/>
   <Atom name="C5" type="1792"/>
   <Atom name="C6" type="1793"/>
   <Atom name="O6" type="1794"/>
   <Atom name="N1" type="1795"/>
   <Atom name="H1" type="1796"/>
   <Atom name="C2" type="1797"/>
   <Atom name="N2" type="1798"/>
   <Atom name="H21" type="1799"/>
   <Atom name="H22" type="1799"/>
   <Atom name="N3" type="1800"/>
   <Atom name="C4" type="1801"/>
   <Atom name="C3'" type="1802"/>
   <Atom name="H3'" type="1803"/>
   <Atom name="C2'" type="1804"/>
   <Atom name="H2'1" type="1805"/>
   <Atom name="O2'" type="1806"/>
   <Atom name="HO'2" type="1807"/>
   <Atom name="O3'" type="1808"/>
   <Bond from="0" to="1"/>
   <Bond from="1" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="5"/>
   <Bond from="5" to="6"/>
   <Bond from="5" to="7"/>
   <Bond from="5" to="25"/>
   <Bond from="7" to="8"/>
   <Bond from="8" to="9"/>
   <Bond from="8" to="10"/>
   <Bond from="8" to="27"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="24"/>
   <Bond from="11" to="12"/>
   <Bond from="11" to="13"/>
   <Bond from="13" to="14"/>
   <Bond from="14" to="15"/>
   <Bond from="14" to="24"/>
   <Bond from="15" to="16"/>
   <Bond from="15" to="17"/>
   <Bond from="17" to="18"/>
   <Bond from="17" to="19"/>
   <Bond from="19" to="20"/>
   <Bond from="19" to="23"/>
   <Bond from="20" to="21"/>
   <Bond from="20" to="22"/>
   <Bond from="23" to="24"/>
   <Bond from="25" to="26"/>
   <Bond from="25" to="27"/>
   <Bond from="25" to="31"/>
   <Bond from="27" to="28"/>
   <Bond from="27" to="29"/>
   <Bond from="29" to="30"/>
   <ExternalBond from="31"/>
  </Residue>
  <Residue name="RGN">
   <Atom name="H5T" type="1809"/>
   <Atom name="O5'" type="1810"/>
   <Atom name="C5'" type="1811"/>
   <Atom name="H5'1" type="1812"/>
   <Atom name="H5'2" type="1812"/>
   <Atom name="C4'" type="1813"/>
   <Atom name="H4'" type="1814"/>
   <Atom name="O4'" type="1815"/>
   <Atom name="C1'" type="1816"/>
   <Atom name="H1'" type="1817"/>
   <Atom name="N9" type="1818"/>
   <Atom name="C8" type="1819"/>
   <Atom name="H8" type="1820"/>
   <Atom name="N7" type="1821"/>
   <Atom name="C5" type="1822"/>
   <Atom name="C6" type="1823"/>
   <Atom name="O6" type="1824"/>
   <Atom name="N1" type="1825"/>
   <Atom name="H1" type="1826"/>
   <Atom name="C2" type="1827"/>
   <Atom name="N2" type="1828"/>
   <Atom name="H21" type="1829"/>
   <Atom name="H22" type="1829"/>
   <Atom name="N3" type="1830"/>
   <Atom name="C4" type="1831"/>
   <Atom name="C3'" type="1832"/>
   <Atom name="H3'" type="1833"/>
   <Atom name="C2'" type="1834"/>
   <Atom name="H2'1" type="1835"/>
   <Atom name="O2'" type="1836"/>
   <Atom name="HO'2" type="1837"/>
   <Atom name="O3'" type="1838"/>
   <Atom name="H3T" type="1839"/>
   <Bond from="0" to="1"/>
   <Bond from="1" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="5"/>
   <Bond from="5" to="6"/>
   <Bond from="5" to="7"/>
   <Bond from="5" to="25"/>
   <Bond from="7" to="8"/>
   <Bond from="8" to="9"/>
   <Bond from="8" to="10"/>
   <Bond from="8" to="27"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="24"/>
   <Bond from="11" to="12"/>
   <Bond from="11" to="13"/>
   <Bond from="13" to="14"/>
   <Bond from="14" to="15"/>
   <Bond from="14" to="24"/>
   <Bond from="15" to="16"/>
   <Bond from="15" to="17"/>
   <Bond from="17" to="18"/>
   <Bond from="17" to="19"/>
   <Bond from="19" to="20"/>
   <Bond from="19" to="23"/>
   <Bond from="20" to="21"/>
   <Bond from="20" to="22"/>
   <Bond from="23" to="24"/>
   <Bond from="25" to="26"/>
   <Bond from="25" to="27"/>
   <Bond from="25" to="31"/>
   <Bond from="27" to="28"/>
   <Bond from="27" to="29"/>
   <Bond from="29" to="30"/>
   <Bond from="31" to="32"/>
  </Residue>
  <Residue name="RU">
   <Atom name="P" type="1840"/>
   <Atom name="O1P" type="1841"/>
   <Atom name="O2P" type="1842"/>
   <Atom name="O5'" type="1843"/>
   <Atom name="C5'" type="1844"/>
   <Atom name="H5'1" type="1845"/>
   <Atom name="H5'2" type="1845"/>
   <Atom name="C4'" type="1846"/>
   <Atom name="H4'" type="1847"/>
   <Atom name="O4'" type="1848"/>
   <Atom name="C1'" type="1849"/>
   <Atom name="H1'" type="1850"/>
   <Atom name="N1" type="1851"/>
   <Atom name="C6" type="1852"/>
   <Atom name="H6" type="1853"/>
   <Atom name="C5" type="1854"/>
   <Atom name="H5" type="1855"/>
   <Atom name="C4" type="1856"/>
   <Atom name="O4" type="1857"/>
   <Atom name="N3" type="1858"/>
   <Atom name="H3" type="1859"/>
   <Atom name="C2" type="1860"/>
   <Atom name="O2" type="1861"/>
   <Atom name="C3'" type="1862"/>
   <Atom name="H3'" type="1863"/>
   <Atom name="C2'" type="1864"/>
   <Atom name="H2'1" type="1865"/>
   <Atom name="O2'" type="1866"/>
   <Atom name="HO'2" type="1867"/>
   <Atom name="O3'" type="1868"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="0" to="3"/>
   <Bond from="3" to="4"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="7" to="9"/>
   <Bond from="7" to="23"/>
   <Bond from="9" to="10"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="12"/>
   <Bond from="10" to="25"/>
   <Bond from="12" to="13"/>
   <Bond from="12" to="21"/>
   <Bond from="13" to="14"/>
   <Bond from="13" to="15"/>
   <Bond from="15" to="16"/>
   <Bond from="15" to="17"/>
   <Bond from="17" to="18"/>
   <Bond from="17" to="19"/>
   <Bond from="19" to="20"/>
   <Bond from="19" to="21"/>
   <Bond from="21" to="22"/>
   <Bond from="23" to="24"/>
   <Bond from="23" to="25"/>
   <Bond from="23" to="29"/>
   <Bond from="25" to="26"/>
   <Bond from="25" to="27"/>
   <Bond from="27" to="28"/>
   <ExternalBond from="0"/>
   <ExternalBond from="29"/>
  </Residue>
  <Residue name="RU3">
   <Atom name="P" type="1869"/>
   <Atom name="O1P" type="1870"/>
   <Atom name="O2P" type="1871"/>
   <Atom name="O5'" type="1872"/>
   <Atom name="C5'" type="1873"/>
   <Atom name="H5'1" type="1874"/>
   <Atom name="H5'2" type="1874"/>
   <Atom name="C4'" type="1875"/>
   <Atom name="H4'" type="1876"/>
   <Atom name="O4'" type="1877"/>
   <Atom name="C1'" type="1878"/>
   <Atom name="H1'" type="1879"/>
   <Atom name="N1" type="1880"/>
   <Atom name="C6" type="1881"/>
   <Atom name="H6" type="1882"/>
   <Atom name="C5" type="1883"/>
   <Atom name="H5" type="1884"/>
   <Atom name="C4" type="1885"/>
   <Atom name="O4" type="1886"/>
   <Atom name="N3" type="1887"/>
   <Atom name="H3" type="1888"/>
   <Atom name="C2" type="1889"/>
   <Atom name="O2" type="1890"/>
   <Atom name="C3'" type="1891"/>
   <Atom name="H3'" type="1892"/>
   <Atom name="C2'" type="1893"/>
   <Atom name="H2'1" type="1894"/>
   <Atom name="O2'" type="1895"/>
   <Atom name="HO'2" type="1896"/>
   <Atom name="O3'" type="1897"/>
   <Atom name="H3T" type="1898"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="0" to="3"/>
   <Bond from="3" to="4"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="7" to="9"/>
   <Bond from="7" to="23"/>
   <Bond from="9" to="10"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="12"/>
   <Bond from="10" to="25"/>
   <Bond from="12" to="13"/>
   <Bond from="12" to="21"/>
   <Bond from="13" to="14"/>
   <Bond from="13" to="15"/>
   <Bond from="15" to="16"/>
   <Bond from="15" to="17"/>
   <Bond from="17" to="18"/>
   <Bond from="17" to="19"/>
   <Bond from="19" to="20"/>
   <Bond from="19" to="21"/>
   <Bond from="21" to="22"/>
   <Bond from="23" to="24"/>
   <Bond from="23" to="25"/>
   <Bond from="23" to="29"/>
   <Bond from="25" to="26"/>
   <Bond from="25" to="27"/>
   <Bond from="27" to="28"/>
   <Bond from="29" to="30"/>
   <ExternalBond from="0"/>
  </Residue>
  <Residue name="RU5">
   <Atom name="H5T" type="1899"/>
   <Atom name="O5'" type="1900"/>
   <Atom name="C5'" type="1901"/>
   <Atom name="H5'1" type="1902"/>
   <Atom name="H5'2" type="1902"/>
   <Atom name="C4'" type="1903"/>
   <Atom name="H4'" type="1904"/>
   <Atom name="O4'" type="1905"/>
   <Atom name="C1'" type="1906"/>
   <Atom name="H1'" type="1907"/>
   <Atom name="N1" type="1908"/>
   <Atom name="C6" type="1909"/>
   <Atom name="H6" type="1910"/>
   <Atom name="C5" type="1911"/>
   <Atom name="H5" type="1912"/>
   <Atom name="C4" type="1913"/>
   <Atom name="O4" type="1914"/>
   <Atom name="N3" type="1915"/>
   <Atom name="H3" type="1916"/>
   <Atom name="C2" type="1917"/>
   <Atom name="O2" type="1918"/>
   <Atom name="C3'" type="1919"/>
   <Atom name="H3'" type="1920"/>
   <Atom name="C2'" type="1921"/>
   <Atom name="H2'1" type="1922"/>
   <Atom name="O2'" type="1923"/>
   <Atom name="HO'2" type="1924"/>
   <Atom name="O3'" type="1925"/>
   <Bond from="0" to="1"/>
   <Bond from="1" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="5"/>
   <Bond from="5" to="6"/>
   <Bond from="5" to="7"/>
   <Bond from="5" to="21"/>
   <Bond from="7" to="8"/>
   <Bond from="8" to="9"/>
   <Bond from="8" to="10"/>
   <Bond from="8" to="23"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="19"/>
   <Bond from="11" to="12"/>
   <Bond from="11" to="13"/>
   <Bond from="13" to="14"/>
   <Bond from="13" to="15"/>
   <Bond from="15" to="16"/>
   <Bond from="15" to="17"/>
   <Bond from="17" to="18"/>
   <Bond from="17" to="19"/>
   <Bond from="19" to="20"/>
   <Bond from="21" to="22"/>
   <Bond from="21" to="23"/>
   <Bond from="21" to="27"/>
   <Bond from="23" to="24"/>
   <Bond from="23" to="25"/>
   <Bond from="25" to="26"/>
   <ExternalBond from="27"/>
  </Residue>
  <Residue name="RUN">
   <Atom name="H5T" type="1926"/>
   <Atom name="O5'" type="1927"/>
   <Atom name="C5'" type="1928"/>
   <Atom name="H5'1" type="1929"/>
   <Atom name="H5'2" type="1929"/>
   <Atom name="C4'" type="1930"/>
   <Atom name="H4'" type="1931"/>
   <Atom name="O4'" type="1932"/>
   <Atom name="C1'" type="1933"/>
   <Atom name="H1'" type="1934"/>
   <Atom name="N1" type="1935"/>
   <Atom name="C6" type="1936"/>
   <Atom name="H6" type="1937"/>
   <Atom name="C5" type="1938"/>
   <Atom name="H5" type="1939"/>
   <Atom name="C4" type="1940"/>
   <Atom name="O4" type="1941"/>
   <Atom name="N3" type="1942"/>
   <Atom name="H3" type="1943"/>
   <Atom name="C2" type="1944"/>
   <Atom name="O2" type="1945"/>
   <Atom name="C3'" type="1946"/>
   <Atom name="H3'" type="1947"/>
   <Atom name="C2'" type="1948"/>
   <Atom name="H2'1" type="1949"/>
   <Atom name="O2'" type="1950"/>
   <Atom name="HO'2" type="1951"/>
   <Atom name="O3'" type="1952"/>
   <Atom name="H3T" type="1953"/>
   <Bond from="0" to="1"/>
   <Bond from="1" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="5"/>
   <Bond from="5" to="6"/>
   <Bond from="5" to="7"/>
   <Bond from="5" to="21"/>
   <Bond from="7" to="8"/>
   <Bond from="8" to="9"/>
   <Bond from="8" to="10"/>
   <Bond from="8" to="23"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="19"/>
   <Bond from="11" to="12"/>
   <Bond from="11" to="13"/>
   <Bond from="13" to="14"/>
   <Bond from="13" to="15"/>
   <Bond from="15" to="16"/>
   <Bond from="15" to="17"/>
   <Bond from="17" to="18"/>
   <Bond from="17" to="19"/>
   <Bond from="19" to="20"/>
   <Bond from="21" to="22"/>
   <Bond from="21" to="23"/>
   <Bond from="21" to="27"/>
   <Bond from="23" to="24"/>
   <Bond from="23" to="25"/>
   <Bond from="25" to="26"/>
   <Bond from="27" to="28"/>
  </Residue>
  <Residue name="Rb+">
   <Atom name="Rb+" type="1960"/>
  </Residue>
  <Residue name="SER">
   <Atom name="N" type="289"/>
   <Atom name="H" type="290"/>
   <Atom name="CA" type="291"/>
   <Atom name="HA" type="292"/>
   <Atom name="CB" type="293"/>
   <Atom name="HB2" type="294"/>
   <Atom name="HB3" type="294"/>
   <Atom name="OG" type="295"/>
   <Atom name="HG" type="296"/>
   <Atom name="C" type="297"/>
   <Atom name="O" type="298"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="9"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="9" to="10"/>
   <ExternalBond from="0"/>
   <ExternalBond from="9"/>
  </Residue>
  <Residue name="THR">
   <Atom name="N" type="299"/>
   <Atom name="H" type="300"/>
   <Atom name="CA" type="301"/>
   <Atom name="HA" type="302"/>
   <Atom name="CB" type="303"/>
   <Atom name="HB" type="304"/>
   <Atom name="CG2" type="305"/>
   <Atom name="HG21" type="306"/>
   <Atom name="HG22" type="306"/>
   <Atom name="HG23" type="306"/>
   <Atom name="OG1" type="307"/>
   <Atom name="HG1" type="308"/>
   <Atom name="C" type="309"/>
   <Atom name="O" type="310"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="12"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="10"/>
   <Bond from="6" to="7"/>
   <Bond from="6" to="8"/>
   <Bond from="6" to="9"/>
   <Bond from="10" to="11"/>
   <Bond from="12" to="13"/>
   <ExternalBond from="0"/>
   <ExternalBond from="12"/>
  </Residue>
  <Residue name="TRP">
   <Atom name="N" type="311"/>
   <Atom name="H" type="312"/>
   <Atom name="CA" type="313"/>
   <Atom name="HA" type="314"/>
   <Atom name="CB" type="315"/>
   <Atom name="HB2" type="316"/>
   <Atom name="HB3" type="316"/>
   <Atom name="CG" type="317"/>
   <Atom name="CD1" type="318"/>
   <Atom name="HD1" type="319"/>
   <Atom name="NE1" type="320"/>
   <Atom name="HE1" type="321"/>
   <Atom name="CE2" type="322"/>
   <Atom name="CZ2" type="323"/>
   <Atom name="HZ2" type="324"/>
   <Atom name="CH2" type="325"/>
   <Atom name="HH2" type="326"/>
   <Atom name="CZ3" type="327"/>
   <Atom name="HZ3" type="328"/>
   <Atom name="CE3" type="329"/>
   <Atom name="HE3" type="330"/>
   <Atom name="CD2" type="331"/>
   <Atom name="C" type="332"/>
   <Atom name="O" type="333"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="22"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="7" to="21"/>
   <Bond from="8" to="9"/>
   <Bond from="8" to="10"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="12"/>
   <Bond from="12" to="13"/>
   <Bond from="12" to="21"/>
   <Bond from="13" to="14"/>
   <Bond from="13" to="15"/>
   <Bond from="15" to="16"/>
   <Bond from="15" to="17"/>
   <Bond from="17" to="18"/>
   <Bond from="17" to="19"/>
   <Bond from="19" to="20"/>
   <Bond from="19" to="21"/>
   <Bond from="22" to="23"/>
   <ExternalBond from="0"/>
   <ExternalBond from="22"/>
  </Residue>
  <Residue name="TYR">
   <Atom name="N" type="334"/>
   <Atom name="H" type="335"/>
   <Atom name="CA" type="336"/>
   <Atom name="HA" type="337"/>
   <Atom name="CB" type="338"/>
   <Atom name="HB2" type="339"/>
   <Atom name="HB3" type="339"/>
   <Atom name="CG" type="340"/>
   <Atom name="CD1" type="341"/>
   <Atom name="HD1" type="342"/>
   <Atom name="CE1" type="343"/>
   <Atom name="HE1" type="344"/>
   <Atom name="CZ" type="345"/>
   <Atom name="OH" type="346"/>
   <Atom name="HH" type="347"/>
   <Atom name="CE2" type="348"/>
   <Atom name="HE2" type="349"/>
   <Atom name="CD2" type="350"/>
   <Atom name="HD2" type="351"/>
   <Atom name="C" type="352"/>
   <Atom name="O" type="353"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="19"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="7"/>
   <Bond from="7" to="8"/>
   <Bond from="7" to="17"/>
   <Bond from="8" to="9"/>
   <Bond from="8" to="10"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="12"/>
   <Bond from="12" to="13"/>
   <Bond from="12" to="15"/>
   <Bond from="13" to="14"/>
   <Bond from="15" to="16"/>
   <Bond from="15" to="17"/>
   <Bond from="17" to="18"/>
   <Bond from="19" to="20"/>
   <ExternalBond from="0"/>
   <ExternalBond from="19"/>
  </Residue>
  <Residue name="VAL">
   <Atom name="N" type="354"/>
   <Atom name="H" type="355"/>
   <Atom name="CA" type="356"/>
   <Atom name="HA" type="357"/>
   <Atom name="CB" type="358"/>
   <Atom name="HB" type="359"/>
   <Atom name="CG1" type="360"/>
   <Atom name="HG11" type="361"/>
   <Atom name="HG12" type="361"/>
   <Atom name="HG13" type="361"/>
   <Atom name="CG2" type="362"/>
   <Atom name="HG21" type="363"/>
   <Atom name="HG22" type="363"/>
   <Atom name="HG23" type="363"/>
   <Atom name="C" type="364"/>
   <Atom name="O" type="365"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
   <Bond from="2" to="3"/>
   <Bond from="2" to="4"/>
   <Bond from="2" to="14"/>
   <Bond from="4" to="5"/>
   <Bond from="4" to="6"/>
   <Bond from="4" to="10"/>
   <Bond from="6" to="7"/>
   <Bond from="6" to="8"/>
   <Bond from="6" to="9"/>
   <Bond from="10" to="11"/>
   <Bond from="10" to="12"/>
   <Bond from="10" to="13"/>
   <Bond from="14" to="15"/>
   <ExternalBond from="0"/>
   <ExternalBond from="14"/>
  </Residue>
 </Residues>
 <HarmonicBondForce>
  <!--Bonds are entered in the same order as the residue list above -->
  <!-- Connections -->
  <!-- ACE-LEU -->
  <Bond class1="C712" class2="N200" length="0.1353" k="312461.12"/>
  <!-- ACE-LEU END -->
  <!-- NME-LEU -->
  <Bond class1="C212" class2="N706" length="0.1353" k="312461.12"/>
  <!-- NME-LEU END -->
  <!-- Connections END-->
  <!-- ACE -->
  <Bond class1="C711" class2="C712" length="0.1523" k="166657.08800000002"/>
  <Bond class1="C711" class2="H710" length="0.10920000000000002" k="273432.76800000004"/>
  <Bond class1="C712" class2="O713" length="0.122" k="581350.064"/>
  <!-- ACE END-->
  <!-- LEU -->
  <Bond class1="N200" class2="H201" length=".1011" k="387907.00800000003"/>
  <Bond class1="N200" class2="C202" length="0.1449" k="211040.96000000002"/>
  <Bond class1="C202" class2="H203" length="0.10920000000000002" k="273432.76800000004"/>
  <Bond class1="C202" class2="C204" length="0.1537" k="159109.152"/>
  <Bond class1="C204" class2="H205" length="0.10920000000000002" k="273432.76800000004"/>
  <Bond class1="C204" class2="C206" length="0.1537" k="159109.152"/>
  <Bond class1="C206" class2="H207" length="0.10920000000000002" k="273432.76800000004"/>
  <Bond class1="C206" class2="C208" length="0.1537" k="159109.152"/>
  <Bond class1="C208" class2="H209" length="0.10920000000000002" k="273432.76800000004"/>
  <Bond class1="C210" class2="C206" length="0.1537" k="159109.152"/>
  <Bond class1="C210" class2="H211" length="0.10920000000000002" k="273432.76800000004"/>
  <Bond class1="C212" class2="C202" length="0.1523" k="166657.08800000002"/>
  <Bond class1="C212" class2="O213" length="0.122" k="581350.064"/>
  <!-- LEU END-->
  <!-- NME -->
  <Bond class1="N706" class2="H707" length="0.1011" k="387907.00800000003"/>
  <Bond class1="N706" class2="C708" length="0.1449" k="211040.96000000002"/>
  <Bond class1="C708" class2="H709" length="0.10920000000000002" k="273432.76800000004"/>
  <!-- NME END -->
 </HarmonicBondForce>
 <HarmonicAngleForce>
  <!-- Angles are entered in the same order as the residue list above -->
  <!-- Connections -->
  <!-- ACE-LEU -->
  <Angle class1="C711" class2="C712" class3="N200" angle="2.018298747006243" k="624.42016"/>
  <Angle class1="C712" class2="N200" class3="C202" angle="2.1273818252558883" k="678.47744"/>
  <Angle class1="C712" class2="N200" class3="H201" angle="2.0799088696016423" k="284.00991999999997"/>
  <Angle class1="N200" class2="C712" class3="O713" angle="2.1432643214490366" k="506.43136000000004"/>
  <!-- ACE-LEU END -->
  <!-- NME-LEU -->
  <Angle class1="N706" class2="C212" class3="O213" angle="2.1432643214490366" k="506.43136000000004"/>
  <Angle class1="N706" class2="C212" class3="C202" angle="2.018298747006243" k="624.42016"/>
  <Angle class1="C212" class2="N706" class3="C708" angle="2.1273818252558883" k="678.47744"/>
  <Angle class1="C212" class2="N706" class3="H707" angle="2.0799088696016423" k="284.00991999999997"/>
  <!-- NME-LEU END -->
  <!-- ConnectionsEND-->
  <!-- ACE -->
  <Angle class1="C711" class2="C712" class3="O713" angle="2.1223203704251046" k="486.85024000000004"/>
  <Angle class1="H710" class2="C711" class3="H710" angle="1.8949039688902434" k="266.85552"/>
  <Angle class1="C712" class2="C711" class3="H710" angle="1.91445165651258" k="435.72176"/>
  <!-- ACE END -->
  <!-- LEU -->
  <Angle class1="N200" class2="C202" class3="H203" angle="1.9109609980085913" k="459.31952"/>
  <Angle class1="N200" class2="C202" class3="C212" angle="1.8919369091618534" k="954.5377599999999"/>
  <Angle class1="H201" class2="N200" class3="C202" angle="2.057219589325716" k="262.67152"/>
  <Angle class1="C202" class2="C212" class3="O213" angle="2.1223203704251046" k="486.85024000000004"/>
  <Angle class1="C202" class2="C204" class3="H205" angle="1.9029324834494175" k="378.98672"/>
  <Angle class1="C202" class2="C204" class3="C206" angle="1.9634954084936207" k="847.00896"/>
  <Angle class1="C204" class2="C206" class3="H207" angle="1.9029324834494175" k="378.98672"/>
  <Angle class1="C204" class2="C206" class3="C210" angle="1.9634954084936207" k="847.00896"/>
  <Angle class1="C204" class2="C206" class3="C208" angle="1.9634954084936207" k="847.00896"/>
  <Angle class1="C204" class2="C202" class3="H203" angle="1.9029324834494175" k="378.98672"/>
  <Angle class1="C204" class2="C202" class3="C212" angle="1.955117828084048" k="961.23216"/>
  <Angle class1="C204" class2="C202" class3="N200" angle="1.9420278586940904" k="795.9641600000001"/>
  <Angle class1="H205" class2="C204" class3="H205" angle="1.8949039688902434" k="266.85552"/>
  <Angle class1="C206" class2="C210" class3="H211" angle="1.9029324834494175" k="378.98672"/>
  <Angle class1="C206" class2="C208" class3="H209" angle="1.9029324834494175" k="378.98672"/>
  <Angle class1="C206" class2="C204" class3="H205" angle="1.9029324834494175" k="378.98672"/>
  <Angle class1="C208" class2="C206" class3="H207" angle="1.9029324834494175" k="378.98672"/>
  <Angle class1="C208" class2="C206" class3="C210" angle="1.9634954084936207" k="847.00896"/>
  <Angle class1="H209" class2="C208" class3="H209" angle="1.8949039688902434" k="266.85552"/>
  <Angle class1="C210" class2="C206" class3="H207" angle="1.9029324834494175" k="378.98672"/>
  <Angle class1="H211" class2="C210" class3="H211" angle="1.8949039688902434" k="266.85552"/>
  <Angle class1="C212" class2="C202" class3="H203" angle="1.91445165651258" k="435.72176"/>
  <!-- LEU END-->
  <!-- NME -->
  <Angle class1="H707" class2="N706" class3="C708" angle="2.057219589325716" k="262.67152"/>
  <Angle class1="H709" class2="C708" class3="N706" angle="1.9109609980085913" k="459.31952"/>
  <Angle class1="H709" class2="C708" class3="H709" angle="1.8949039688902434" k="266.85552"/>
  <!-- NME END-->
 </HarmonicAngleForce>
 <PeriodicTorsionForce>
  <!-- Torsions are in the same residue order as above -->
  <!-- Connections -->
  <!-- ACE-LEU -->
  <Proper class1="C711" class2="C712" class3="N200" class4="H201" k1="0" k2="10.250800000000002" k3="0" k4="0" periodicity1="1" periodicity2="2" periodicity3="3" periodicity4="4" phase1="0" phase2="3.141592653589793" phase3="0" phase4="3.141592653589793"/>
  <Proper class1="C711" class2="C712" class3="N200" class4="C202" k1="4.811599999999999" k2="12.738188000000001" k3="0" k4="0" periodicity1="1" periodicity2="2" periodicity3="3" periodicity4="4" phase1="0" phase2="3.141592653589793" phase3="0" phase4="3.141592653589793"/>
  <Proper class1="H710" class2="C711" class3="C712" class4="N200" k1="0" k2="0" k3="0" k4="0" periodicity1="1" periodicity2="2" periodicity3="3" periodicity4="4" phase1="0" phase2="3.141592653589793" phase3="0" phase4="3.141592653589793"/>
  <Proper class1="H710" class2="C711" class3="C712" class4="O713" k1="0" k2="0" k3="0" k4="0" periodicity1="1" periodicity2="2" periodicity3="3" periodicity4="4" phase1="0" phase2="3.141592653589793" phase3="0" phase4="3.141592653589793"/>
  <Proper class1="C712" class2="N200" class3="C202" class4="C204" k1="-.7740400000000001" k2="-0.33472" k3="0" k4="0" periodicity1="1" periodicity2="2" periodicity3="3" periodicity4="4" phase1="0" phase2="3.141592653589793" phase3="0" phase4="3.141592653589793"/>
  <Proper class1="C712" class2="N200" class3="C202" class4="C212" k1=".04184" k2="0.8577199999999999" k3="-0.18828" k4="0" periodicity1="1" periodicity2="2" periodicity3="3" periodicity4="4" phase1="0" phase2="3.141592653589793" phase3="0" phase4="3.141592653589793"/>
  <Proper class1="C712" class2="N200" class3="C202" class4="H203" k1="0" k2="0" k3="0" k4="0" periodicity1="1" periodicity2="2" periodicity3="3" periodicity4="4" phase1="0" phase2="3.141592653589793" phase3="0" phase4="3.141592653589793"/>
  <Proper class1="O713" class2="C712" class3="N200" class4="H201" k1="0" k2="10.250800000000002" k3="0" k4="0" periodicity1="1" periodicity2="2" periodicity3="3" periodicity4="4" phase1="0" phase2="3.141592653589793" phase3="0" phase4="3.141592653589793"/>
  <Proper class1="O713" class2="C712" class3="N200" class4="C202" k1="0" k2="12.738188000000001" k3="0" k4="0" periodicity1="1" periodicity2="2" periodicity3="3" periodicity4="4" phase1="0" phase2="3.141592653589793" phase3="0" phase4="3.141592653589793"/>
  <Proper class1="H201" class2="N200" class3="C712" class4="C202" k1="0" k2="10.46" k3="0" k4="0" periodicity1="1" periodicity2="2" periodicity3="3" periodicity4="4" phase1="0" phase2="3.141592653589793" phase3="0" phase4="3.141592653589793"/>
  <Improper class1="C712" class2="C711" class3="N200" class4="O713" k1="0" k2="43.932" k3="0" k4="0" periodicity1="1" periodicity2="2" periodicity3="3" periodicity4="4" phase1="0" phase2="3.141592653589793" phase3="0" phase4="3.141592653589793"/>
  <Improper class1="N200" class2="C202" class3="C712" class4="H201" k1="0" k2="10.46" k3="0" k4="0" periodicity1="1" periodicity2="2" periodicity3="3" periodicity4="4" phase1="0" phase2="3.141592653589793" phase3="0" phase4="3.141592653589793"/>

  <!-- ACE-LEU END -->
  <!-- NME-LEU -->
  <Proper class1="N706" class2="C212" class3="C202" class4="C204" k1="-.5439200000000001" k2="-0.66944" k3="-0.50208" k4="0" periodicity1="1" periodicity2="2" periodicity3="3" periodicity4="4" phase1="0" phase2="3.141592653589793" phase3="0" phase4="3.141592653589793"/>
  <Proper class1="N706" class2="C212" class3="C202" class4="N200" k1=".18828" k2="0.89956" k3="-0.60668" k4="0" periodicity1="1" periodicity2="2" periodicity3="3" periodicity4="4" phase1="0" phase2="3.141592653589793" phase3="0" phase4="3.141592653589793"/>
  <Proper class1="N706" class2="C212" class3="C202" class4="H203" k1="0" k2="0" k3="0" k4="0" periodicity1="1" periodicity2="2" periodicity3="3" periodicity4="4" phase1="0" phase2="3.141592653589793" phase3="0" phase4="3.141592653589793"/>
  <Proper class1="H707" class2="N706" class3="C212" class4="C202" k1="0" k2="10.250800000000002" k3="0" k4="0" periodicity1="1" periodicity2="2" periodicity3="3" periodicity4="4" phase1="0" phase2="3.141592653589793" phase3="0" phase4="3.141592653589793"/>
  <Proper class1="H707" class2="N706" class3="C212" class4="O213" k1="0" k2="10.250800000000002" k3="0" k4="0" periodicity1="1" periodicity2="2" periodicity3="3" periodicity4="4" phase1="0" phase2="3.141592653589793" phase3="0" phase4="3.141592653589793"/>
  <Proper class1="C708" class2="N706" class3="C212" class4="C202" k1="4.811599999999999" k2="12.738188000000001" k3="0" k4="0" periodicity1="1" periodicity2="2" periodicity3="3" periodicity4="4" phase1="0" phase2="3.141592653589793" phase3="0" phase4="3.141592653589793"/>
  <Proper class1="C708" class2="N706" class3="C212" class4="O213" k1="0" k2="12.738188000000001" k3="0" k4="0" periodicity1="1" periodicity2="2" periodicity3="3" periodicity4="4" phase1="0" phase2="3.141592653589793" phase3="0" phase4="3.141592653589793"/>
  <Proper class1="H709" class2="C708" class3="N706" class4="C212" k1="0" k2="0" k3="0" k4="0" periodicity1="1" periodicity2="2" periodicity3="3" periodicity4="4" phase1="0" phase2="3.141592653589793" phase3="0" phase4="3.141592653589793"/>
  <Improper class1="C212" class2="C202" class3="N706" class4="O213" k1="0" k2="43.932" k3="0" k4="0" periodicity1="1" periodicity2="2" periodicity3="3" periodicity4="4" phase1="0" phase2="3.141592653589793" phase3="0" phase4="3.141592653589793"/>
  <Improper class1="N706" class2="C708" class3="C212" class4="H707" k1="0" k2="10.46" k3="0" k4="0" periodicity1="1" periodicity2="2" periodicity3="3" periodicity4="4" phase1="0" phase2="3.141592653589793" phase3="0" phase4="3.141592653589793"/>
  <Proper class1="H707" class2="N706" class3="C212" class4="C708" k1="0" k2="10.46" k3="0" k4="0" periodicity1="1" periodicity2="2" periodicity3="3" periodicity4="4" phase1="0" phase2="3.141592653589793" phase3="0" phase4="3.141592653589793"/>
  <!-- NME-LEU END -->
  <!-- Connections END -->
  <!-- ACE -->
  <Proper class1="H710" class2="C711" class3="C712" class4="O713" periodicity1="1" phase1="0" k1="0" periodicity2="2" phase2="0" k2="0" periodicity3="3" phase3="0" k3="0" periodicity4="4" phase4="3.141592653589793" k4="0"/>
  <!-- ACE END-->
  <!-- LEU -->
  <Proper class1="N200" class2="C202" class3="C204" class4="C206" k1="-4.351360000000001" k2="-1.3598000000000001" k3="1.21336" k4="0" periodicity1="1" periodicity2="2" periodicity3="3" periodicity4="4" phase1="0" phase2="3.141592653589793" phase3="0" phase4="3.141592653589793"/>
  <Proper class1="N200" class2="C202" class3="C204" class4="H205" k1="0" k2="0" k3=".9706880000000001" k4="0" periodicity1="1" periodicity2="2" periodicity3="3" periodicity4="4" phase1="0" phase2="3.141592653589793" phase3="0" phase4="3.141592653589793"/>
  <Proper class1="N200" class2="C202" class3="C212" class4="O213" k1="0" k2="0" k3="0" k4="0" periodicity1="1" periodicity2="2" periodicity3="3" periodicity4="4" phase1="0" phase2="3.141592653589793" phase3="0" phase4="3.141592653589793"/>
  <Proper class1="H201" class2="N200" class3="C202" class4="C204" k1="0" k2="0" k3="0" k4="0" periodicity1="1" periodicity2="2" periodicity3="3" periodicity4="4" phase1="0" phase2="3.141592653589793" phase3="0" phase4="3.141592653589793"/>
  <Proper class1="H201" class2="N200" class3="C202" class4="C212" k1="0" k2="0" k3="0" k4="0" periodicity1="1" periodicity2="2" periodicity3="3" periodicity4="4" phase1="0" phase2="3.141592653589793" phase3="0" phase4="3.141592653589793"/>
  <Proper class1="H201" class2="N200" class3="C202" class4="H203" k1="0" k2="0" k3="0" k4="0" periodicity1="1" periodicity2="2" periodicity3="3" periodicity4="4" phase1="0" phase2="3.141592653589793" phase3="0" phase4="3.141592653589793"/>
  <Proper class1="C202" class2="C204" class3="C206" class4="C208" k1="-2.2802800000000003" k2="1.48532" k3=".41840000000000005" k4="0" periodicity1="1" periodicity2="2" periodicity3="3" periodicity4="4" phase1="0" phase2="3.141592653589793" phase3="0" phase4="3.141592653589793"/>
  <Proper class1="C202" class2="C204" class3="C206" class4="C210" k1="-2.2802800000000003" k2="1.48532" k3=".41840000000000005" k4="0" periodicity1="1" periodicity2="2" periodicity3="3" periodicity4="4" phase1="0" phase2="3.141592653589793" phase3="0" phase4="3.141592653589793"/>
  <Proper class1="C202" class2="C204" class3="C206" class4="H207" k1="0" k2="0" k3=".6276" k4="0" periodicity1="1" periodicity2="2" periodicity3="3" periodicity4="4" phase1="0" phase2="3.141592653589793" phase3="0" phase4="3.141592653589793"/>
  <Proper class1="H203" class2="C202" class3="C204" class4="C206" k1="0" k2="0" k3=".6276" k4="0" periodicity1="1" periodicity2="2" periodicity3="3" periodicity4="4" phase1="0" phase2="3.141592653589793" phase3="0" phase4="3.141592653589793"/>
  <Proper class1="H203" class2="C202" class3="C204" class4="H205" k1="0" k2="0" k3=".6276" k4="0" periodicity1="1" periodicity2="2" periodicity3="3" periodicity4="4" phase1="0" phase2="3.141592653589793" phase3="0" phase4="3.141592653589793"/>
  <Proper class1="H203" class2="C202" class3="C212" class4="O213" k1="0" k2="0" k3="0" k4="0" periodicity1="1" periodicity2="2" periodicity3="3" periodicity4="4" phase1="0" phase2="3.141592653589793" phase3="0" phase4="3.141592653589793"/>
  <Proper class1="C204" class2="C206" class3="C208" class4="H209" k1="0" k2="0" k3=".6276" k4="0" periodicity1="1" periodicity2="2" periodicity3="3" periodicity4="4" phase1="0" phase2="3.141592653589793" phase3="0" phase4="3.141592653589793"/>
  <Proper class1="C204" class2="C206" class3="C210" class4="H211" k1="0" k2="0" k3=".6276" k4="0" periodicity1="1" periodicity2="2" periodicity3="3" periodicity4="4" phase1="0" phase2="3.141592653589793" phase3="0" phase4="3.141592653589793"/>
  <Proper class1="C204" class2="C202" class3="C212" class4="O213" k1="0" k2="0" k3="0" k4="0" periodicity1="1" periodicity2="2" periodicity3="3" periodicity4="4" phase1="0" phase2="3.141592653589793" phase3="0" phase4="3.141592653589793"/>
  <Proper class1="H205" class2="C204" class3="C202" class4="C212" k1="0" k2="0" k3="-.158992" k4="0" periodicity1="1" periodicity2="2" periodicity3="3" periodicity4="4" phase1="0" phase2="3.141592653589793" phase3="0" phase4="3.141592653589793"/>
  <Proper class1="H205" class2="C204" class3="C206" class4="C210" k1="0" k2="0" k3=".6276" k4="0" periodicity1="1" periodicity2="2" periodicity3="3" periodicity4="4" phase1="0" phase2="3.141592653589793" phase3="0" phase4="3.141592653589793"/>
  <Proper class1="H205" class2="C204" class3="C206" class4="H207" k1="0" k2="0" k3=".6276" k4="0" periodicity1="1" periodicity2="2" periodicity3="3" periodicity4="4" phase1="0" phase2="3.141592653589793" phase3="0" phase4="3.141592653589793"/>
  <Proper class1="H205" class2="C204" class3="C202" class4="C212" k1="0" k2="0" k3="-.158992" k4="0" periodicity1="1" periodicity2="2" periodicity3="3" periodicity4="4" phase1="0" phase2="3.141592653589793" phase3="0" phase4="3.141592653589793"/>
  <Proper class1="H205" class2="C204" class3="C206" class4="C208" k1="0" k2="0" k3=".6276" k4="0" periodicity1="1" periodicity2="2" periodicity3="3" periodicity4="4" phase1="0" phase2="3.141592653589793" phase3="0" phase4="3.141592653589793"/>
  <Proper class1="C206" class2="C204" class3="C202" class4="C212" k1="1.1087600000000002" k2="0.7112800000000001" k3="0.25104" k4="0" periodicity1="1" periodicity2="2" periodicity3="3" periodicity4="4" phase1="0" phase2="3.141592653589793" phase3="0" phase4="3.141592653589793"/>
  <Proper class1="H207" class2="C206" class3="C208" class4="H209" k1="0" k2="0" k3=".6276" k4="0" periodicity1="1" periodicity2="2" periodicity3="3" periodicity4="4" phase1="0" phase2="3.141592653589793" phase3="0" phase4="3.141592653589793"/>
  <Proper class1="H207" class2="C206" class3="C210" class4="H211" k1="0" k2="0" k3=".6276" k4="0" periodicity1="1" periodicity2="2" periodicity3="3" periodicity4="4" phase1="0" phase2="3.141592653589793" phase3="0" phase4="3.141592653589793"/>
  <Proper class1="C208" class2="C206" class3="C210" class4="H211" k1="0" k2="0" k3=".6276" k4="0" periodicity1="1" periodicity2="2" periodicity3="3" periodicity4="4" phase1="0" phase2="3.141592653589793" phase3="0" phase4="3.141592653589793"/>
  <Proper class1="H209" class2="C208" class3="C206" class4="C210" k1="0" k2="0" k3=".6276" k4="0" periodicity1="1" periodicity2="2" periodicity3="3" periodicity4="4" phase1="0" phase2="3.141592653589793" phase3="0" phase4="3.141592653589793"/>
  <!-- LEU END -->
  <!-- NME -->
  <Proper class1="H707" class2="N706" class3="C708" class4="H709" periodicity1="1" phase1="0" k1="0" periodicity2="2" phase2="0" k2="0" periodicity3="3" phase3="0" k3="0" periodicity4="4" phase4="3.141592653589793" k4="0"/>
 <!-- NME END-->
 </PeriodicTorsionForce>
 <NonbondedForce coulomb14scale="0.5" lj14scale="0.5">
  <Atom type="0" charge="-0.4157" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1" charge="0.2719" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="2" charge="0.0337" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="3" charge="0.0823" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="4" charge="-0.1825" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="5" charge="0.0603" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="6" charge="0.5973" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="7" charge="-0.5679" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="8" charge="-0.3479" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="9" charge="0.2747" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="10" charge="-0.2637" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="11" charge="0.156" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="12" charge="-0.0007" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="13" charge="0.0327" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="14" charge="0.039" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="15" charge="0.0285" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="16" charge="0.0486" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="17" charge="0.0687" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="18" charge="-0.5295" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="19" charge="0.3456" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="20" charge="0.8076" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="21" charge="-0.8627" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="22" charge="0.4478" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="23" charge="-0.8627" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="24" charge="0.4478" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="25" charge="0.7341" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="26" charge="-0.5894" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="27" charge="-0.4157" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="28" charge="0.2719" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="29" charge="0.0341" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="30" charge="0.0864" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="31" charge="-0.0316" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="32" charge="0.0488" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="33" charge="0.6462" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="34" charge="-0.5554" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="35" charge="-0.6376" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="36" charge="0.4747" sigma="0.0" epsilon="0.0"/>
  <Atom type="37" charge="0.5973" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="38" charge="-0.5679" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="39" charge="-0.4157" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="40" charge="0.2719" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="41" charge="0.0143" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="42" charge="0.1048" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="43" charge="-0.2041" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="44" charge="0.0797" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="45" charge="0.713" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="46" charge="-0.5931" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="47" charge="-0.9191" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="48" charge="0.4196" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="49" charge="0.5973" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="50" charge="-0.5679" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="51" charge="-0.5163" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="52" charge="0.2936" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="53" charge="0.0381" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="54" charge="0.088" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="55" charge="-0.0303" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="56" charge="-0.0122" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="57" charge="0.7994" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="58" charge="-0.8014" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="59" charge="-0.8014" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="60" charge="0.5366" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="61" charge="-0.5819" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="62" charge="-0.4157" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="63" charge="0.2719" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="64" charge="-0.0351" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="65" charge="0.0508" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="66" charge="-0.2413" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="67" charge="0.1122" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="68" charge="-0.8844" sigma="0.356359487256" epsilon="1.046"/>
  <Atom type="69" charge="0.5973" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="70" charge="-0.5679" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="71" charge="-0.4157" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="72" charge="0.2719" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="73" charge="0.0213" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="74" charge="0.1124" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="75" charge="-0.1231" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="76" charge="0.1112" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="77" charge="-0.3119" sigma="0.356359487256" epsilon="1.046"/>
  <Atom type="78" charge="0.1933" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="79" charge="0.5973" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="80" charge="-0.5679" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="81" charge="-0.4157" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="82" charge="0.2719" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="83" charge="0.0429" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="84" charge="0.0766" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="85" charge="-0.079" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="86" charge="0.091" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="87" charge="-0.1081" sigma="0.356359487256" epsilon="1.046"/>
  <Atom type="88" charge="0.5973" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="89" charge="-0.5679" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="90" charge="-0.4157" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="91" charge="0.2719" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="92" charge="0.0145" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="93" charge="0.0779" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="94" charge="-0.0071" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="95" charge="0.0256" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="96" charge="-0.0174" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="97" charge="0.043" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="98" charge="0.6801" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="99" charge="-0.5838" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="100" charge="-0.6511" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="101" charge="0.4641" sigma="0.0" epsilon="0.0"/>
  <Atom type="102" charge="0.5973" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="103" charge="-0.5679" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="104" charge="-0.4157" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="105" charge="0.2719" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="106" charge="-0.0031" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="107" charge="0.085" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="108" charge="-0.0036" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="109" charge="0.0171" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="110" charge="-0.0645" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="111" charge="0.0352" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="112" charge="0.6951" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="113" charge="-0.6086" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="114" charge="-0.9407" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="115" charge="0.4251" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="116" charge="0.5973" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="117" charge="-0.5679" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="118" charge="-0.5163" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="119" charge="0.2936" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="120" charge="0.0397" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="121" charge="0.1105" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="122" charge="0.056" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="123" charge="-0.0173" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="124" charge="0.0136" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="125" charge="-0.0425" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="126" charge="0.8054" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="127" charge="-0.8188" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="128" charge="-0.8188" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="129" charge="0.5366" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="130" charge="-0.5819" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="131" charge="-0.4157" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="132" charge="0.2719" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="133" charge="-0.0252" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="134" charge="0.0698" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="135" charge="0.5973" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="136" charge="-0.5679" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="137" charge="-0.4157" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="138" charge="0.2719" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="139" charge="0.0188" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="140" charge="0.0881" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="141" charge="-0.0462" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="142" charge="0.0402" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="143" charge="-0.0266" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="144" charge="-0.3811" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="145" charge="0.3649" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="146" charge="0.2057" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="147" charge="0.1392" sigma="0.242146271591" epsilon="0.06276"/>
  <Atom type="148" charge="-0.5727" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="149" charge="0.1292" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="150" charge="0.1147" sigma="0.251055258772" epsilon="0.06276"/>
  <Atom type="151" charge="0.5973" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="152" charge="-0.5679" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="153" charge="-0.4157" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="154" charge="0.2719" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="155" charge="-0.0581" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="156" charge="0.136" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="157" charge="-0.0074" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="158" charge="0.0367" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="159" charge="0.1868" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="160" charge="-0.5432" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="161" charge="0.1635" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="162" charge="0.1435" sigma="0.242146271591" epsilon="0.06276"/>
  <Atom type="163" charge="-0.2795" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="164" charge="0.3339" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="165" charge="-0.2207" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="166" charge="0.1862" sigma="0.251055258772" epsilon="0.06276"/>
  <Atom type="167" charge="0.5973" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="168" charge="-0.5679" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="169" charge="-0.3479" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="170" charge="0.2747" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="171" charge="-0.1354" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="172" charge="0.1212" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="173" charge="-0.0414" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="174" charge="0.081" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="175" charge="-0.0012" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="176" charge="-0.1513" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="177" charge="0.3866" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="178" charge="-0.017" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="179" charge="0.2681" sigma="0.242146271591" epsilon="0.06276"/>
  <Atom type="180" charge="-0.1718" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="181" charge="0.3911" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="182" charge="-0.1141" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="183" charge="0.2317" sigma="0.251055258772" epsilon="0.06276"/>
  <Atom type="184" charge="0.7341" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="185" charge="-0.5894" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="186" charge="-0.4157" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="187" charge="0.2719" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="188" charge="-0.0597" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="189" charge="0.0869" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="190" charge="0.1303" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="191" charge="0.0187" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="192" charge="-0.3204" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="193" charge="0.0882" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="194" charge="-0.043" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="195" charge="0.0236" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="196" charge="-0.066" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="197" charge="0.0186" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="198" charge="0.5973" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="199" charge="-0.5679" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="200" charge="1" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="201" charge="1" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="202" charge="1" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="203" charge="1" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="204" charge="1" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="205" charge="1" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="206" charge="1" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="207" charge="1" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="208" charge="1" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="209" charge="1" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="210" charge="1" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="211" charge="1" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="212" charge="1" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="213" charge="1" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="214" charge="1" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="215" charge="1" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="216" charge="1" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="217" charge="1" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="218" charge="1" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="219" charge="1" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="220" charge="1" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="221" charge="1" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="222" charge="1" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="223" charge="1" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="224" charge="1" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="225" charge="1" sigma="0.195997717991" epsilon="0.0656888"/>
  <Atom type="226" charge="1" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="227" charge="1" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="228" charge="1" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="229" charge="1" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="230" charge="1" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="231" charge="1" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="232" charge="1" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="233" charge="0.1426" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="234" charge="-0.0094" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="235" charge="0.0362" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="236" charge="0.0187" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="237" charge="0.0103" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="238" charge="-0.0479" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="239" charge="0.0621" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="240" charge="-0.0143" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="241" charge="0.1135" sigma="0.195997717991" epsilon="0.0656888"/>
  <Atom type="242" charge="-0.3854" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="243" charge="0.34" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="244" charge="0.7341" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="245" charge="-0.5894" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="246" charge="-0.4157" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="247" charge="0.2719" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="248" charge="-0.0237" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="249" charge="0.088" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="250" charge="0.0342" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="251" charge="0.0241" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="252" charge="0.0018" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="253" charge="0.044" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="254" charge="-0.2737" sigma="0.356359487256" epsilon="1.046"/>
  <Atom type="255" charge="-0.0536" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="256" charge="0.0684" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="257" charge="0.5973" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="258" charge="-0.5679" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="259" charge="-0.4157" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="260" charge="0.2719" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="261" charge="-0.0024" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="262" charge="0.0978" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="263" charge="-0.0343" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="264" charge="0.0295" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="265" charge="0.0118" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="266" charge="-0.1256" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="267" charge="0.133" sigma="0.259964245953" epsilon="0.06276"/>
  <Atom type="268" charge="-0.1704" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="269" charge="0.143" sigma="0.259964245953" epsilon="0.06276"/>
  <Atom type="270" charge="-0.1072" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="271" charge="0.1297" sigma="0.259964245953" epsilon="0.06276"/>
  <Atom type="272" charge="-0.1704" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="273" charge="0.143" sigma="0.259964245953" epsilon="0.06276"/>
  <Atom type="274" charge="-0.1256" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="275" charge="0.133" sigma="0.259964245953" epsilon="0.06276"/>
  <Atom type="276" charge="0.5973" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="277" charge="-0.5679" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="278" charge="-0.2548" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="279" charge="0.0192" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="280" charge="0.0391" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="281" charge="0.0189" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="282" charge="0.0213" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="283" charge="-0.007" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="284" charge="0.0253" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="285" charge="-0.0266" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="286" charge="0.0641" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="287" charge="0.5896" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="288" charge="-0.5748" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="289" charge="-0.4157" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="290" charge="0.2719" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="291" charge="-0.0249" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="292" charge="0.0843" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="293" charge="0.2117" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="294" charge="0.0352" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="295" charge="-0.6546" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="296" charge="0.4275" sigma="0.0" epsilon="0.0"/>
  <Atom type="297" charge="0.5973" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="298" charge="-0.5679" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="299" charge="-0.4157" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="300" charge="0.2719" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="301" charge="-0.0389" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="302" charge="0.1007" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="303" charge="0.3654" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="304" charge="0.0043" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="305" charge="-0.2438" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="306" charge="0.0642" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="307" charge="-0.6761" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="308" charge="0.4102" sigma="0.0" epsilon="0.0"/>
  <Atom type="309" charge="0.5973" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="310" charge="-0.5679" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="311" charge="-0.4157" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="312" charge="0.2719" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="313" charge="-0.0275" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="314" charge="0.1123" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="315" charge="-0.005" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="316" charge="0.0339" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="317" charge="-0.1415" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="318" charge="-0.1638" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="319" charge="0.2062" sigma="0.251055258772" epsilon="0.06276"/>
  <Atom type="320" charge="-0.3418" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="321" charge="0.3412" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="322" charge="0.138" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="323" charge="-0.2601" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="324" charge="0.1572" sigma="0.259964245953" epsilon="0.06276"/>
  <Atom type="325" charge="-0.1134" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="326" charge="0.1417" sigma="0.259964245953" epsilon="0.06276"/>
  <Atom type="327" charge="-0.1972" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="328" charge="0.1447" sigma="0.259964245953" epsilon="0.06276"/>
  <Atom type="329" charge="-0.2387" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="330" charge="0.17" sigma="0.259964245953" epsilon="0.06276"/>
  <Atom type="331" charge="0.1243" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="332" charge="0.5973" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="333" charge="-0.5679" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="334" charge="-0.4157" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="335" charge="0.2719" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="336" charge="-0.0014" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="337" charge="0.0876" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="338" charge="-0.0152" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="339" charge="0.0295" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="340" charge="-0.0011" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="341" charge="-0.1906" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="342" charge="0.1699" sigma="0.259964245953" epsilon="0.06276"/>
  <Atom type="343" charge="-0.2341" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="344" charge="0.1656" sigma="0.259964245953" epsilon="0.06276"/>
  <Atom type="345" charge="0.3226" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="346" charge="-0.5579" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="347" charge="0.3992" sigma="0.0" epsilon="0.0"/>
  <Atom type="348" charge="-0.2341" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="349" charge="0.1656" sigma="0.259964245953" epsilon="0.06276"/>
  <Atom type="350" charge="-0.1906" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="351" charge="0.1699" sigma="0.259964245953" epsilon="0.06276"/>
  <Atom type="352" charge="0.5973" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="353" charge="-0.5679" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="354" charge="-0.4157" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="355" charge="0.2719" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="356" charge="-0.0875" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="357" charge="0.0969" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="358" charge="0.2985" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="359" charge="-0.0297" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="360" charge="-0.3192" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="361" charge="0.0791" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="362" charge="-0.3192" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="363" charge="0.0791" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="364" charge="0.5973" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="365" charge="-0.5679" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="366" charge="-0.3821" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="367" charge="0.2681" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="368" charge="-0.1747" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="369" charge="0.1067" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="370" charge="-0.2093" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="371" charge="0.0764" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="372" charge="0.7731" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="373" charge="-0.8055" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="374" charge="-0.8055" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="375" charge="-0.3481" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="376" charge="0.2764" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="377" charge="-0.3068" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="378" charge="0.1447" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="379" charge="-0.0374" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="380" charge="0.0371" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="381" charge="0.0744" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="382" charge="0.0185" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="383" charge="0.1114" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="384" charge="0.0468" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="385" charge="-0.5564" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="386" charge="0.3479" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="387" charge="0.8368" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="388" charge="-0.8737" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="389" charge="0.4493" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="390" charge="-0.8737" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="391" charge="0.4493" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="392" charge="0.8557" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="393" charge="-0.8266" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="394" charge="-0.8266" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="395" charge="-0.3821" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="396" charge="0.2681" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="397" charge="-0.208" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="398" charge="0.1358" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="399" charge="-0.2299" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="400" charge="0.1023" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="401" charge="0.7153" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="402" charge="-0.601" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="403" charge="-0.9084" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="404" charge="0.415" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="405" charge="0.805" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="406" charge="-0.8147" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="407" charge="-0.8147" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="408" charge="-0.5192" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="409" charge="0.3055" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="410" charge="-0.1817" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="411" charge="0.1046" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="412" charge="-0.0677" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="413" charge="-0.0212" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="414" charge="0.8851" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="415" charge="-0.8162" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="416" charge="-0.8162" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="417" charge="0.7256" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="418" charge="-0.7887" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="419" charge="-0.7887" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="420" charge="-0.3821" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="421" charge="0.2681" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="422" charge="-0.1635" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="423" charge="0.1396" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="424" charge="-0.1996" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="425" charge="0.1437" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="426" charge="-0.3102" sigma="0.356359487256" epsilon="1.046"/>
  <Atom type="427" charge="0.2068" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="428" charge="0.7497" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="429" charge="-0.7981" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="430" charge="-0.7981" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="431" charge="-0.3821" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="432" charge="0.2681" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="433" charge="-0.1318" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="434" charge="0.0938" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="435" charge="-0.1943" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="436" charge="0.1228" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="437" charge="-0.0529" sigma="0.356359487256" epsilon="1.046"/>
  <Atom type="438" charge="0.7618" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="439" charge="-0.8041" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="440" charge="-0.8041" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="441" charge="-0.3821" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="442" charge="0.2681" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="443" charge="-0.2248" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="444" charge="0.1232" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="445" charge="-0.0664" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="446" charge="0.0452" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="447" charge="-0.021" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="448" charge="0.0203" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="449" charge="0.7093" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="450" charge="-0.6098" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="451" charge="-0.9574" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="452" charge="0.4304" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="453" charge="0.7775" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="454" charge="-0.8042" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="455" charge="-0.8042" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="456" charge="-0.5192" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="457" charge="0.3055" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="458" charge="-0.2059" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="459" charge="0.1399" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="460" charge="0.0071" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="461" charge="-0.0078" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="462" charge="0.0675" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="463" charge="-0.0548" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="464" charge="0.8183" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="465" charge="-0.822" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="466" charge="-0.822" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="467" charge="0.742" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="468" charge="-0.793" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="469" charge="-0.793" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="470" charge="-0.3821" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="471" charge="0.2681" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="472" charge="-0.2493" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="473" charge="0.1056" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="474" charge="0.7231" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="475" charge="-0.7855" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="476" charge="-0.7855" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="477" charge="-0.3821" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="478" charge="0.2681" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="479" charge="-0.1739" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="480" charge="0.11" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="481" charge="-0.1046" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="482" charge="0.0565" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="483" charge="0.0293" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="484" charge="-0.3892" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="485" charge="0.3755" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="486" charge="0.1925" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="487" charge="0.1418" sigma="0.242146271591" epsilon="0.06276"/>
  <Atom type="488" charge="-0.5629" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="489" charge="0.1001" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="490" charge="0.1241" sigma="0.251055258772" epsilon="0.06276"/>
  <Atom type="491" charge="0.7615" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="492" charge="-0.8016" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="493" charge="-0.8016" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="494" charge="-0.3821" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="495" charge="0.2681" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="496" charge="-0.2699" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="497" charge="0.165" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="498" charge="-0.1068" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="499" charge="0.062" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="500" charge="0.2724" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="501" charge="-0.5517" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="502" charge="0.1558" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="503" charge="0.1448" sigma="0.242146271591" epsilon="0.06276"/>
  <Atom type="504" charge="-0.267" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="505" charge="0.3319" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="506" charge="-0.2588" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="507" charge="0.1957" sigma="0.251055258772" epsilon="0.06276"/>
  <Atom type="508" charge="0.7916" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="509" charge="-0.8065" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="510" charge="-0.8065" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="511" charge="-0.3481" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="512" charge="0.2764" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="513" charge="-0.1445" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="514" charge="0.1115" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="515" charge="-0.08" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="516" charge="0.0868" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="517" charge="0.0298" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="518" charge="-0.1501" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="519" charge="0.3883" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="520" charge="-0.0251" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="521" charge="0.2694" sigma="0.242146271591" epsilon="0.06276"/>
  <Atom type="522" charge="-0.1683" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="523" charge="0.3913" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="524" charge="-0.1256" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="525" charge="0.2336" sigma="0.251055258772" epsilon="0.06276"/>
  <Atom type="526" charge="0.8032" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="527" charge="-0.8177" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="528" charge="-0.8177" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="529" charge="-0.3821" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="530" charge="0.2681" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="531" charge="-0.31" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="532" charge="0.1375" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="533" charge="0.0363" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="534" charge="0.0766" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="535" charge="-0.3498" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="536" charge="0.1021" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="537" charge="-0.0323" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="538" charge="0.0321" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="539" charge="-0.0699" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="540" charge="0.0196" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="541" charge="0.8343" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="542" charge="-0.819" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="543" charge="-0.819" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="544" charge="-0.3821" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="545" charge="0.2681" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="546" charge="-0.2847" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="547" charge="0.1346" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="548" charge="-0.2469" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="549" charge="0.0974" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="550" charge="0.3706" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="551" charge="-0.0374" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="552" charge="-0.4163" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="553" charge="0.1038" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="554" charge="-0.4163" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="555" charge="0.1038" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="556" charge="0.8326" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="557" charge="-0.8199" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="558" charge="-0.8199" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="559" charge="-0.3481" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="560" charge="0.2764" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="561" charge="-0.2903" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="562" charge="0.1438" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="563" charge="-0.0538" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="564" charge="0.0482" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="565" charge="0.0227" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="566" charge="0.0134" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="567" charge="-0.0392" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="568" charge="0.0611" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="569" charge="-0.0176" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="570" charge="0.1121" sigma="0.195997717991" epsilon="0.0656888"/>
  <Atom type="571" charge="-0.3741" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="572" charge="0.3374" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="573" charge="0.8488" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="574" charge="-0.8252" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="575" charge="-0.8252" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="576" charge="-0.3821" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="577" charge="0.2681" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="578" charge="-0.2597" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="579" charge="0.1277" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="580" charge="-0.0236" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="581" charge="0.048" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="582" charge="0.0492" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="583" charge="0.0317" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="584" charge="-0.2692" sigma="0.356359487256" epsilon="1.046"/>
  <Atom type="585" charge="-0.0376" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="586" charge="0.0625" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="587" charge="0.8013" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="588" charge="-0.8105" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="589" charge="-0.8105" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="590" charge="-0.3821" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="591" charge="0.2681" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="592" charge="-0.1825" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="593" charge="0.1098" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="594" charge="-0.0959" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="595" charge="0.0443" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="596" charge="0.0552" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="597" charge="-0.13" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="598" charge="0.1408" sigma="0.259964245953" epsilon="0.06276"/>
  <Atom type="599" charge="-0.1847" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="600" charge="0.1461" sigma="0.259964245953" epsilon="0.06276"/>
  <Atom type="601" charge="-0.0944" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="602" charge="0.128" sigma="0.259964245953" epsilon="0.06276"/>
  <Atom type="603" charge="-0.1847" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="604" charge="0.1461" sigma="0.259964245953" epsilon="0.06276"/>
  <Atom type="605" charge="-0.13" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="606" charge="0.1408" sigma="0.259964245953" epsilon="0.06276"/>
  <Atom type="607" charge="0.766" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="608" charge="-0.8026" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="609" charge="-0.8026" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="610" charge="-0.2802" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="611" charge="0.0434" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="612" charge="0.0331" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="613" charge="0.0466" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="614" charge="0.0172" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="615" charge="-0.0543" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="616" charge="0.0381" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="617" charge="-0.1336" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="618" charge="0.0776" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="619" charge="0.6631" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="620" charge="-0.7697" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="621" charge="-0.7697" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="622" charge="-0.3821" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="623" charge="0.2681" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="624" charge="-0.2722" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="625" charge="0.1304" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="626" charge="0.1123" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="627" charge="0.0813" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="628" charge="-0.6514" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="629" charge="0.4474" sigma="0.0" epsilon="0.0"/>
  <Atom type="630" charge="0.8113" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="631" charge="-0.8132" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="632" charge="-0.8132" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="633" charge="-0.3821" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="634" charge="0.2681" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="635" charge="-0.242" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="636" charge="0.1207" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="637" charge="0.3025" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="638" charge="0.0078" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="639" charge="-0.1853" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="640" charge="0.0586" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="641" charge="-0.6496" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="642" charge="0.4119" sigma="0.0" epsilon="0.0"/>
  <Atom type="643" charge="0.781" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="644" charge="-0.8044" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="645" charge="-0.8044" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="646" charge="-0.3821" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="647" charge="0.2681" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="648" charge="-0.2084" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="649" charge="0.1272" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="650" charge="-0.0742" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="651" charge="0.0497" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="652" charge="-0.0796" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="653" charge="-0.1808" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="654" charge="0.2043" sigma="0.251055258772" epsilon="0.06276"/>
  <Atom type="655" charge="-0.3316" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="656" charge="0.3413" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="657" charge="0.1222" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="658" charge="-0.2594" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="659" charge="0.1567" sigma="0.259964245953" epsilon="0.06276"/>
  <Atom type="660" charge="-0.102" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="661" charge="0.1401" sigma="0.259964245953" epsilon="0.06276"/>
  <Atom type="662" charge="-0.2287" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="663" charge="0.1507" sigma="0.259964245953" epsilon="0.06276"/>
  <Atom type="664" charge="-0.1837" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="665" charge="0.1491" sigma="0.259964245953" epsilon="0.06276"/>
  <Atom type="666" charge="0.1078" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="667" charge="0.7658" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="668" charge="-0.8011" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="669" charge="-0.8011" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="670" charge="-0.3821" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="671" charge="0.2681" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="672" charge="-0.2015" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="673" charge="0.1092" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="674" charge="-0.0752" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="675" charge="0.049" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="676" charge="0.0243" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="677" charge="-0.1922" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="678" charge="0.178" sigma="0.259964245953" epsilon="0.06276"/>
  <Atom type="679" charge="-0.2458" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="680" charge="0.1673" sigma="0.259964245953" epsilon="0.06276"/>
  <Atom type="681" charge="0.3395" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="682" charge="-0.5643" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="683" charge="0.4017" sigma="0.0" epsilon="0.0"/>
  <Atom type="684" charge="-0.2458" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="685" charge="0.1673" sigma="0.259964245953" epsilon="0.06276"/>
  <Atom type="686" charge="-0.1922" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="687" charge="0.178" sigma="0.259964245953" epsilon="0.06276"/>
  <Atom type="688" charge="0.7817" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="689" charge="-0.807" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="690" charge="-0.807" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="691" charge="-0.3821" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="692" charge="0.2681" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="693" charge="-0.3438" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="694" charge="0.1438" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="695" charge="0.194" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="696" charge="0.0308" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="697" charge="-0.3064" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="698" charge="0.0836" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="699" charge="1" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="700" charge="1" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="701" charge="1" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="702" charge="1" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="703" charge="1" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="704" charge="1" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="705" charge="1" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="706" charge="1" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="707" charge="1" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="708" charge="1" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="709" charge="1" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="710" charge="1" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="711" charge="1" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="712" charge="1" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="713" charge="1" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="714" charge="1" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="715" charge="1" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="716" charge="1" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="717" charge="1" sigma="0.195997717991" epsilon="0.0656888"/>
  <Atom type="718" charge="-0.0597" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="719" charge="0.03" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="720" charge="0.6163" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="721" charge="-0.5722" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="722" charge="0.1305" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="723" charge="0.2083" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="724" charge="-0.0223" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="725" charge="0.1242" sigma="0.195997717991" epsilon="0.0656888"/>
  <Atom type="726" charge="0.0118" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="727" charge="0.0226" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="728" charge="0.0236" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="729" charge="0.0309" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="730" charge="0.0935" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="731" charge="0.0527" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="732" charge="-0.565" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="733" charge="0.3592" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="734" charge="0.8281" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="735" charge="-0.8693" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="736" charge="0.4494" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="737" charge="-0.8693" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="738" charge="0.4494" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="739" charge="0.7214" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="740" charge="-0.6013" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="741" charge="0.1801" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="742" charge="0.1921" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="743" charge="0.0368" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="744" charge="0.1231" sigma="0.195997717991" epsilon="0.0656888"/>
  <Atom type="745" charge="-0.0283" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="746" charge="0.0515" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="747" charge="0.5833" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="748" charge="-0.5744" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="749" charge="-0.8634" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="750" charge="0.4097" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="751" charge="0.6163" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="752" charge="-0.5722" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="753" charge="0.0782" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="754" charge="0.22" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="755" charge="0.0292" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="756" charge="0.1141" sigma="0.195997717991" epsilon="0.0656888"/>
  <Atom type="757" charge="-0.0235" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="758" charge="-0.0169" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="759" charge="0.8194" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="760" charge="-0.8084" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="761" charge="-0.8084" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="762" charge="0.5621" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="763" charge="-0.5889" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="764" charge="0.1325" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="765" charge="0.2023" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="766" charge="0.0927" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="767" charge="0.1411" sigma="0.195997717991" epsilon="0.0656888"/>
  <Atom type="768" charge="-0.1195" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="769" charge="0.1188" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="770" charge="-0.3298" sigma="0.356359487256" epsilon="1.046"/>
  <Atom type="771" charge="0.1975" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="772" charge="0.6123" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="773" charge="-0.5713" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="774" charge="0.2069" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="775" charge="0.1815" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="776" charge="0.1055" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="777" charge="0.0922" sigma="0.195997717991" epsilon="0.0656888"/>
  <Atom type="778" charge="-0.0277" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="779" charge="0.068" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="780" charge="-0.0984" sigma="0.356359487256" epsilon="1.046"/>
  <Atom type="781" charge="0.6123" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="782" charge="-0.5713" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="783" charge="0.1493" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="784" charge="0.1996" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="785" charge="0.0536" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="786" charge="0.1015" sigma="0.195997717991" epsilon="0.0656888"/>
  <Atom type="787" charge="0.0651" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="788" charge="0.005" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="789" charge="-0.0903" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="790" charge="0.0331" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="791" charge="0.7354" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="792" charge="-0.6133" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="793" charge="-1.0031" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="794" charge="0.4429" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="795" charge="0.6123" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="796" charge="-0.5713" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="797" charge="0.0017" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="798" charge="0.2391" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="799" charge="0.0588" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="800" charge="0.1202" sigma="0.195997717991" epsilon="0.0656888"/>
  <Atom type="801" charge="0.0909" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="802" charge="-0.0232" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="803" charge="-0.0236" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="804" charge="-0.0315" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="805" charge="0.8087" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="806" charge="-0.8189" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="807" charge="-0.8189" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="808" charge="0.5621" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="809" charge="-0.5889" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="810" charge="0.2943" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="811" charge="0.1642" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="812" charge="-0.01" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="813" charge="0.0895" sigma="0.195997717991" epsilon="0.0656888"/>
  <Atom type="814" charge="0.6163" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="815" charge="-0.5722" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="816" charge="0.1542" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="817" charge="0.1963" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="818" charge="0.0964" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="819" charge="0.0958" sigma="0.195997717991" epsilon="0.0656888"/>
  <Atom type="820" charge="0.0259" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="821" charge="0.0209" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="822" charge="-0.0399" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="823" charge="-0.3819" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="824" charge="0.3632" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="825" charge="0.2127" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="826" charge="0.1385" sigma="0.242146271591" epsilon="0.06276"/>
  <Atom type="827" charge="-0.5711" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="828" charge="0.1046" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="829" charge="0.1299" sigma="0.251055258772" epsilon="0.06276"/>
  <Atom type="830" charge="0.6123" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="831" charge="-0.5713" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="832" charge="0.1472" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="833" charge="0.2016" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="834" charge="0.0236" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="835" charge="0.138" sigma="0.195997717991" epsilon="0.0656888"/>
  <Atom type="836" charge="0.0489" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="837" charge="0.0223" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="838" charge="0.174" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="839" charge="-0.5579" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="840" charge="0.1804" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="841" charge="0.1397" sigma="0.242146271591" epsilon="0.06276"/>
  <Atom type="842" charge="-0.2781" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="843" charge="0.3324" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="844" charge="-0.2349" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="845" charge="0.1963" sigma="0.251055258772" epsilon="0.06276"/>
  <Atom type="846" charge="0.6123" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="847" charge="-0.5713" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="848" charge="0.256" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="849" charge="0.1704" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="850" charge="0.0581" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="851" charge="0.1047" sigma="0.195997717991" epsilon="0.0656888"/>
  <Atom type="852" charge="0.0484" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="853" charge="0.0531" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="854" charge="-0.0236" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="855" charge="-0.151" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="856" charge="0.3821" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="857" charge="-0.0011" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="858" charge="0.2645" sigma="0.242146271591" epsilon="0.06276"/>
  <Atom type="859" charge="-0.1739" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="860" charge="0.3921" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="861" charge="-0.1433" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="862" charge="0.2495" sigma="0.251055258772" epsilon="0.06276"/>
  <Atom type="863" charge="0.7214" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="864" charge="-0.6013" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="865" charge="0.0311" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="866" charge="0.2329" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="867" charge="0.0257" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="868" charge="0.1031" sigma="0.195997717991" epsilon="0.0656888"/>
  <Atom type="869" charge="0.1885" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="870" charge="0.0213" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="871" charge="-0.372" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="872" charge="0.0947" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="873" charge="-0.0387" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="874" charge="0.0201" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="875" charge="-0.0908" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="876" charge="0.0226" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="877" charge="0.6123" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="878" charge="-0.5713" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="879" charge="0.101" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="880" charge="0.2148" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="881" charge="0.0104" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="882" charge="0.1053" sigma="0.195997717991" epsilon="0.0656888"/>
  <Atom type="883" charge="-0.0244" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="884" charge="0.0256" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="885" charge="0.3421" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="886" charge="-0.038" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="887" charge="-0.4106" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="888" charge="0.098" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="889" charge="-0.4104" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="890" charge="0.098" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="891" charge="0.6123" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="892" charge="-0.5713" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="893" charge="0.0966" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="894" charge="0.2165" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="895" charge="-0.0015" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="896" charge="0.118" sigma="0.195997717991" epsilon="0.0656888"/>
  <Atom type="897" charge="0.0212" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="898" charge="0.0283" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="899" charge="-0.0048" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="900" charge="0.0121" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="901" charge="-0.0608" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="902" charge="0.0633" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="903" charge="-0.0181" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="904" charge="0.1171" sigma="0.195997717991" epsilon="0.0656888"/>
  <Atom type="905" charge="-0.3764" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="906" charge="0.3382" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="907" charge="0.7214" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="908" charge="-0.6013" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="909" charge="0.1592" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="910" charge="0.1984" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="911" charge="0.0221" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="912" charge="0.1116" sigma="0.195997717991" epsilon="0.0656888"/>
  <Atom type="913" charge="0.0865" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="914" charge="0.0125" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="915" charge="0.0334" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="916" charge="0.0292" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="917" charge="-0.2774" sigma="0.356359487256" epsilon="1.046"/>
  <Atom type="918" charge="-0.0341" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="919" charge="0.0597" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="920" charge="0.6123" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="921" charge="-0.5713" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="922" charge="0.1737" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="923" charge="0.1921" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="924" charge="0.0733" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="925" charge="0.1041" sigma="0.195997717991" epsilon="0.0656888"/>
  <Atom type="926" charge="0.033" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="927" charge="0.0104" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="928" charge="0.0031" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="929" charge="-0.1392" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="930" charge="0.1374" sigma="0.259964245953" epsilon="0.06276"/>
  <Atom type="931" charge="-0.1602" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="932" charge="0.1433" sigma="0.259964245953" epsilon="0.06276"/>
  <Atom type="933" charge="-0.1208" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="934" charge="0.1329" sigma="0.259964245953" epsilon="0.06276"/>
  <Atom type="935" charge="-0.1603" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="936" charge="0.1433" sigma="0.259964245953" epsilon="0.06276"/>
  <Atom type="937" charge="-0.1391" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="938" charge="0.1374" sigma="0.259964245953" epsilon="0.06276"/>
  <Atom type="939" charge="0.6123" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="940" charge="-0.5713" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="941" charge="-0.202" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="942" charge="0.312" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="943" charge="-0.012" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="944" charge="0.1" sigma="0.195997717991" epsilon="0.0656888"/>
  <Atom type="945" charge="-0.121" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="946" charge="0.1" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="947" charge="-0.115" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="948" charge="0.1" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="949" charge="0.1" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="950" charge="0.1" sigma="0.195997717991" epsilon="0.0656888"/>
  <Atom type="951" charge="0.526" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="952" charge="-0.5" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="953" charge="0.1849" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="954" charge="0.1898" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="955" charge="0.0567" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="956" charge="0.0782" sigma="0.195997717991" epsilon="0.0656888"/>
  <Atom type="957" charge="0.2596" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="958" charge="0.0273" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="959" charge="-0.6714" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="960" charge="0.4239" sigma="0.0" epsilon="0.0"/>
  <Atom type="961" charge="0.6163" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="962" charge="-0.5722" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="963" charge="0.1812" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="964" charge="0.1934" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="965" charge="0.0034" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="966" charge="0.1087" sigma="0.195997717991" epsilon="0.0656888"/>
  <Atom type="967" charge="0.4514" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="968" charge="-0.0323" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="969" charge="-0.2554" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="970" charge="0.0627" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="971" charge="-0.6764" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="972" charge="0.407" sigma="0.0" epsilon="0.0"/>
  <Atom type="973" charge="0.6163" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="974" charge="-0.5722" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="975" charge="0.1913" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="976" charge="0.1888" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="977" charge="0.0421" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="978" charge="0.1162" sigma="0.195997717991" epsilon="0.0656888"/>
  <Atom type="979" charge="0.0543" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="980" charge="0.0222" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="981" charge="-0.1654" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="982" charge="-0.1788" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="983" charge="0.2195" sigma="0.251055258772" epsilon="0.06276"/>
  <Atom type="984" charge="-0.3444" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="985" charge="0.3412" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="986" charge="0.1575" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="987" charge="-0.271" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="988" charge="0.1589" sigma="0.259964245953" epsilon="0.06276"/>
  <Atom type="989" charge="-0.108" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="990" charge="0.1411" sigma="0.259964245953" epsilon="0.06276"/>
  <Atom type="991" charge="-0.2034" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="992" charge="0.1458" sigma="0.259964245953" epsilon="0.06276"/>
  <Atom type="993" charge="-0.2265" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="994" charge="0.1646" sigma="0.259964245953" epsilon="0.06276"/>
  <Atom type="995" charge="0.1132" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="996" charge="0.6123" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="997" charge="-0.5713" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="998" charge="0.194" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="999" charge="0.1873" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="1000" charge="0.057" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1001" charge="0.0983" sigma="0.195997717991" epsilon="0.0656888"/>
  <Atom type="1002" charge="0.0659" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1003" charge="0.0102" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="1004" charge="-0.0205" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1005" charge="-0.2002" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1006" charge="0.172" sigma="0.259964245953" epsilon="0.06276"/>
  <Atom type="1007" charge="-0.2239" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1008" charge="0.165" sigma="0.259964245953" epsilon="0.06276"/>
  <Atom type="1009" charge="0.3139" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1010" charge="-0.5578" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="1011" charge="0.4001" sigma="0.0" epsilon="0.0"/>
  <Atom type="1012" charge="-0.2239" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1013" charge="0.165" sigma="0.259964245953" epsilon="0.06276"/>
  <Atom type="1014" charge="-0.2002" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1015" charge="0.172" sigma="0.259964245953" epsilon="0.06276"/>
  <Atom type="1016" charge="0.6123" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1017" charge="-0.5713" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1018" charge="0.0577" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1019" charge="0.2272" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="1020" charge="-0.0054" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1021" charge="0.1093" sigma="0.195997717991" epsilon="0.0656888"/>
  <Atom type="1022" charge="0.3196" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1023" charge="-0.0221" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="1024" charge="-0.3129" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1025" charge="0.0735" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="1026" charge="-0.3129" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1027" charge="0.0735" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="1028" charge="0.6163" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1029" charge="-0.5722" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1030" charge="1.1659" sigma="0.374177461619" epsilon="0.8368"/>
  <Atom type="1031" charge="-0.7761" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1032" charge="-0.7761" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1033" charge="-0.4954" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1034" charge="-0.0069" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1035" charge="0.0754" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1036" charge="0.1629" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1037" charge="0.1176" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1038" charge="-0.3691" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1039" charge="0.0431" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1040" charge="0.1838" sigma="0.229317330049" epsilon="0.0656888"/>
  <Atom type="1041" charge="-0.0268" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1042" charge="0.1607" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1043" charge="0.1877" sigma="0.242146271591" epsilon="0.06276"/>
  <Atom type="1044" charge="-0.6175" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1045" charge="0.0725" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1046" charge="0.6897" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1047" charge="-0.9123" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1048" charge="0.4167" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="1049" charge="-0.7624" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1050" charge="0.5716" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1051" charge="0.0598" sigma="0.242146271591" epsilon="0.06276"/>
  <Atom type="1052" charge="-0.7417" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1053" charge="0.38" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1054" charge="0.0713" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1055" charge="0.0985" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1056" charge="-0.0854" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1057" charge="0.0718" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="1058" charge="-0.5232" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1059" charge="1.1659" sigma="0.374177461619" epsilon="0.8368"/>
  <Atom type="1060" charge="-0.7761" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1061" charge="-0.7761" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1062" charge="-0.4954" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1063" charge="-0.0069" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1064" charge="0.0754" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1065" charge="0.1629" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1066" charge="0.1176" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1067" charge="-0.3691" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1068" charge="0.0431" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1069" charge="0.1838" sigma="0.229317330049" epsilon="0.0656888"/>
  <Atom type="1070" charge="-0.0268" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1071" charge="0.1607" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1072" charge="0.1877" sigma="0.242146271591" epsilon="0.06276"/>
  <Atom type="1073" charge="-0.6175" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1074" charge="0.0725" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1075" charge="0.6897" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1076" charge="-0.9123" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1077" charge="0.4167" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="1078" charge="-0.7624" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1079" charge="0.5716" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1080" charge="0.0598" sigma="0.242146271591" epsilon="0.06276"/>
  <Atom type="1081" charge="-0.7417" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1082" charge="0.38" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1083" charge="0.0713" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1084" charge="0.0985" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1085" charge="-0.0854" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1086" charge="0.0718" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="1087" charge="-0.6549" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="1088" charge="0.4396" sigma="0.0" epsilon="0.0"/>
  <Atom type="1089" charge="0.4422" sigma="0.0" epsilon="0.0"/>
  <Atom type="1090" charge="-0.6318" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="1091" charge="-0.0069" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1092" charge="0.0754" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1093" charge="0.1629" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1094" charge="0.1176" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1095" charge="-0.3691" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1096" charge="0.0431" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1097" charge="0.1838" sigma="0.229317330049" epsilon="0.0656888"/>
  <Atom type="1098" charge="-0.0268" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1099" charge="0.1607" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1100" charge="0.1877" sigma="0.242146271591" epsilon="0.06276"/>
  <Atom type="1101" charge="-0.6175" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1102" charge="0.0725" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1103" charge="0.6897" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1104" charge="-0.9123" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1105" charge="0.4167" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="1106" charge="-0.7624" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1107" charge="0.5716" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1108" charge="0.0598" sigma="0.242146271591" epsilon="0.06276"/>
  <Atom type="1109" charge="-0.7417" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1110" charge="0.38" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1111" charge="0.0713" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1112" charge="0.0985" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1113" charge="-0.0854" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1114" charge="0.0718" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="1115" charge="-0.5232" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1116" charge="0.4422" sigma="0.0" epsilon="0.0"/>
  <Atom type="1117" charge="-0.6318" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="1118" charge="-0.0069" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1119" charge="0.0754" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1120" charge="0.1629" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1121" charge="0.1176" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1122" charge="-0.3691" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1123" charge="0.0431" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1124" charge="0.1838" sigma="0.229317330049" epsilon="0.0656888"/>
  <Atom type="1125" charge="-0.0268" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1126" charge="0.1607" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1127" charge="0.1877" sigma="0.242146271591" epsilon="0.06276"/>
  <Atom type="1128" charge="-0.6175" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1129" charge="0.0725" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1130" charge="0.6897" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1131" charge="-0.9123" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1132" charge="0.4167" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="1133" charge="-0.7624" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1134" charge="0.5716" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1135" charge="0.0598" sigma="0.242146271591" epsilon="0.06276"/>
  <Atom type="1136" charge="-0.7417" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1137" charge="0.38" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1138" charge="0.0713" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1139" charge="0.0985" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1140" charge="-0.0854" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1141" charge="0.0718" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="1142" charge="-0.6549" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="1143" charge="0.4396" sigma="0.0" epsilon="0.0"/>
  <Atom type="1144" charge="1.1659" sigma="0.374177461619" epsilon="0.8368"/>
  <Atom type="1145" charge="-0.7761" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1146" charge="-0.7761" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1147" charge="-0.4954" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1148" charge="-0.0069" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1149" charge="0.0754" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1150" charge="0.1629" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1151" charge="0.1176" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1152" charge="-0.3691" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1153" charge="-0.0116" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1154" charge="0.1963" sigma="0.229317330049" epsilon="0.0656888"/>
  <Atom type="1155" charge="-0.0339" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1156" charge="-0.0183" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1157" charge="0.2293" sigma="0.251055258772" epsilon="0.06276"/>
  <Atom type="1158" charge="-0.5222" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1159" charge="0.1863" sigma="0.259964245953" epsilon="0.06276"/>
  <Atom type="1160" charge="0.8439" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1161" charge="-0.9773" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1162" charge="0.4314" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="1163" charge="-0.7748" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1164" charge="0.7959" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1165" charge="-0.6548" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1166" charge="0.0713" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1167" charge="0.0985" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1168" charge="-0.0854" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1169" charge="0.0718" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="1170" charge="-0.5232" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1171" charge="1.1659" sigma="0.374177461619" epsilon="0.8368"/>
  <Atom type="1172" charge="-0.7761" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1173" charge="-0.7761" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1174" charge="-0.4954" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1175" charge="-0.0069" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1176" charge="0.0754" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1177" charge="0.1629" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1178" charge="0.1176" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1179" charge="-0.3691" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1180" charge="-0.0116" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1181" charge="0.1963" sigma="0.229317330049" epsilon="0.0656888"/>
  <Atom type="1182" charge="-0.0339" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1183" charge="-0.0183" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1184" charge="0.2293" sigma="0.251055258772" epsilon="0.06276"/>
  <Atom type="1185" charge="-0.5222" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1186" charge="0.1863" sigma="0.259964245953" epsilon="0.06276"/>
  <Atom type="1187" charge="0.8439" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1188" charge="-0.9773" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1189" charge="0.4314" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="1190" charge="-0.7748" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1191" charge="0.7959" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1192" charge="-0.6548" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1193" charge="0.0713" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1194" charge="0.0985" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1195" charge="-0.0854" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1196" charge="0.0718" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="1197" charge="-0.6549" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="1198" charge="0.4396" sigma="0.0" epsilon="0.0"/>
  <Atom type="1199" charge="0.4422" sigma="0.0" epsilon="0.0"/>
  <Atom type="1200" charge="-0.6318" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="1201" charge="-0.0069" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1202" charge="0.0754" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1203" charge="0.1629" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1204" charge="0.1176" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1205" charge="-0.3691" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1206" charge="-0.0116" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1207" charge="0.1963" sigma="0.229317330049" epsilon="0.0656888"/>
  <Atom type="1208" charge="-0.0339" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1209" charge="-0.0183" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1210" charge="0.2293" sigma="0.251055258772" epsilon="0.06276"/>
  <Atom type="1211" charge="-0.5222" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1212" charge="0.1863" sigma="0.259964245953" epsilon="0.06276"/>
  <Atom type="1213" charge="0.8439" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1214" charge="-0.9773" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1215" charge="0.4314" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="1216" charge="-0.7748" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1217" charge="0.7959" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1218" charge="-0.6548" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1219" charge="0.0713" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1220" charge="0.0985" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1221" charge="-0.0854" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1222" charge="0.0718" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="1223" charge="-0.5232" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1224" charge="0.4422" sigma="0.0" epsilon="0.0"/>
  <Atom type="1225" charge="-0.6318" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="1226" charge="-0.0069" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1227" charge="0.0754" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1228" charge="0.1629" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1229" charge="0.1176" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1230" charge="-0.3691" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1231" charge="-0.0116" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1232" charge="0.1963" sigma="0.229317330049" epsilon="0.0656888"/>
  <Atom type="1233" charge="-0.0339" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1234" charge="-0.0183" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1235" charge="0.2293" sigma="0.251055258772" epsilon="0.06276"/>
  <Atom type="1236" charge="-0.5222" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1237" charge="0.1863" sigma="0.259964245953" epsilon="0.06276"/>
  <Atom type="1238" charge="0.8439" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1239" charge="-0.9773" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1240" charge="0.4314" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="1241" charge="-0.7748" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1242" charge="0.7959" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1243" charge="-0.6548" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1244" charge="0.0713" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1245" charge="0.0985" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1246" charge="-0.0854" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1247" charge="0.0718" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="1248" charge="-0.6549" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="1249" charge="0.4396" sigma="0.0" epsilon="0.0"/>
  <Atom type="1250" charge="1.1659" sigma="0.374177461619" epsilon="0.8368"/>
  <Atom type="1251" charge="-0.7761" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1252" charge="-0.7761" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1253" charge="-0.4954" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1254" charge="-0.0069" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1255" charge="0.0754" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1256" charge="0.1629" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1257" charge="0.1176" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1258" charge="-0.3691" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1259" charge="0.0358" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1260" charge="0.1746" sigma="0.229317330049" epsilon="0.0656888"/>
  <Atom type="1261" charge="0.0577" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1262" charge="0.0736" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1263" charge="0.1997" sigma="0.242146271591" epsilon="0.06276"/>
  <Atom type="1264" charge="-0.5725" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1265" charge="0.1991" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1266" charge="0.4918" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1267" charge="-0.5699" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1268" charge="-0.5053" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1269" charge="0.352" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="1270" charge="0.7432" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1271" charge="-0.923" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1272" charge="0.4235" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="1273" charge="-0.6636" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1274" charge="0.1814" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1275" charge="0.0713" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1276" charge="0.0985" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1277" charge="-0.0854" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1278" charge="0.0718" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="1279" charge="-0.5232" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1280" charge="1.1659" sigma="0.374177461619" epsilon="0.8368"/>
  <Atom type="1281" charge="-0.7761" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1282" charge="-0.7761" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1283" charge="-0.4954" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1284" charge="-0.0069" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1285" charge="0.0754" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1286" charge="0.1629" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1287" charge="0.1176" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1288" charge="-0.3691" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1289" charge="0.0358" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1290" charge="0.1746" sigma="0.229317330049" epsilon="0.0656888"/>
  <Atom type="1291" charge="0.0577" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1292" charge="0.0736" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1293" charge="0.1997" sigma="0.242146271591" epsilon="0.06276"/>
  <Atom type="1294" charge="-0.5725" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1295" charge="0.1991" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1296" charge="0.4918" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1297" charge="-0.5699" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1298" charge="-0.5053" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1299" charge="0.352" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="1300" charge="0.7432" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1301" charge="-0.923" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1302" charge="0.4235" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="1303" charge="-0.6636" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1304" charge="0.1814" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1305" charge="0.0713" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1306" charge="0.0985" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1307" charge="-0.0854" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1308" charge="0.0718" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="1309" charge="-0.6549" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="1310" charge="0.4396" sigma="0.0" epsilon="0.0"/>
  <Atom type="1311" charge="0.4422" sigma="0.0" epsilon="0.0"/>
  <Atom type="1312" charge="-0.6318" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="1313" charge="-0.0069" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1314" charge="0.0754" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1315" charge="0.1629" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1316" charge="0.1176" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1317" charge="-0.3691" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1318" charge="0.0358" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1319" charge="0.1746" sigma="0.229317330049" epsilon="0.0656888"/>
  <Atom type="1320" charge="0.0577" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1321" charge="0.0736" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1322" charge="0.1997" sigma="0.242146271591" epsilon="0.06276"/>
  <Atom type="1323" charge="-0.5725" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1324" charge="0.1991" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1325" charge="0.4918" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1326" charge="-0.5699" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1327" charge="-0.5053" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1328" charge="0.352" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="1329" charge="0.7432" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1330" charge="-0.923" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1331" charge="0.4235" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="1332" charge="-0.6636" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1333" charge="0.1814" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1334" charge="0.0713" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1335" charge="0.0985" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1336" charge="-0.0854" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1337" charge="0.0718" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="1338" charge="-0.5232" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1339" charge="0.4422" sigma="0.0" epsilon="0.0"/>
  <Atom type="1340" charge="-0.6318" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="1341" charge="-0.0069" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1342" charge="0.0754" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1343" charge="0.1629" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1344" charge="0.1176" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1345" charge="-0.3691" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1346" charge="0.0358" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1347" charge="0.1746" sigma="0.229317330049" epsilon="0.0656888"/>
  <Atom type="1348" charge="0.0577" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1349" charge="0.0736" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1350" charge="0.1997" sigma="0.242146271591" epsilon="0.06276"/>
  <Atom type="1351" charge="-0.5725" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1352" charge="0.1991" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1353" charge="0.4918" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1354" charge="-0.5699" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1355" charge="-0.5053" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1356" charge="0.352" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="1357" charge="0.7432" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1358" charge="-0.923" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1359" charge="0.4235" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="1360" charge="-0.6636" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1361" charge="0.1814" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1362" charge="0.0713" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1363" charge="0.0985" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1364" charge="-0.0854" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1365" charge="0.0718" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="1366" charge="-0.6549" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="1367" charge="0.4396" sigma="0.0" epsilon="0.0"/>
  <Atom type="1368" charge="1.1659" sigma="0.374177461619" epsilon="0.8368"/>
  <Atom type="1369" charge="-0.7761" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1370" charge="-0.7761" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1371" charge="-0.4954" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1372" charge="-0.0069" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1373" charge="0.0754" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1374" charge="0.1629" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1375" charge="0.1176" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1376" charge="-0.3691" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1377" charge="0.068" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1378" charge="0.1804" sigma="0.229317330049" epsilon="0.0656888"/>
  <Atom type="1379" charge="-0.0239" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1380" charge="-0.2209" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1381" charge="0.2607" sigma="0.251055258772" epsilon="0.06276"/>
  <Atom type="1382" charge="0.0025" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1383" charge="-0.2269" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1384" charge="0.077" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="1385" charge="0.5194" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1386" charge="-0.5563" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1387" charge="-0.434" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1388" charge="0.342" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="1389" charge="0.5677" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1390" charge="-0.5881" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1391" charge="0.0713" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1392" charge="0.0985" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1393" charge="-0.0854" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1394" charge="0.0718" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="1395" charge="-0.5232" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1396" charge="1.1659" sigma="0.374177461619" epsilon="0.8368"/>
  <Atom type="1397" charge="-0.7761" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1398" charge="-0.7761" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1399" charge="-0.4954" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1400" charge="-0.0069" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1401" charge="0.0754" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1402" charge="0.1629" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1403" charge="0.1176" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1404" charge="-0.3691" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1405" charge="0.068" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1406" charge="0.1804" sigma="0.229317330049" epsilon="0.0656888"/>
  <Atom type="1407" charge="-0.0239" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1408" charge="-0.2209" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1409" charge="0.2607" sigma="0.251055258772" epsilon="0.06276"/>
  <Atom type="1410" charge="0.0025" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1411" charge="-0.2269" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1412" charge="0.077" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="1413" charge="0.5194" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1414" charge="-0.5563" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1415" charge="-0.434" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1416" charge="0.342" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="1417" charge="0.5677" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1418" charge="-0.5881" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1419" charge="0.0713" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1420" charge="0.0985" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1421" charge="-0.0854" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1422" charge="0.0718" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="1423" charge="-0.6549" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="1424" charge="0.4396" sigma="0.0" epsilon="0.0"/>
  <Atom type="1425" charge="0.4422" sigma="0.0" epsilon="0.0"/>
  <Atom type="1426" charge="-0.6318" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="1427" charge="-0.0069" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1428" charge="0.0754" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1429" charge="0.1629" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1430" charge="0.1176" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1431" charge="-0.3691" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1432" charge="0.068" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1433" charge="0.1804" sigma="0.229317330049" epsilon="0.0656888"/>
  <Atom type="1434" charge="-0.0239" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1435" charge="-0.2209" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1436" charge="0.2607" sigma="0.251055258772" epsilon="0.06276"/>
  <Atom type="1437" charge="0.0025" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1438" charge="-0.2269" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1439" charge="0.077" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="1440" charge="0.5194" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1441" charge="-0.5563" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1442" charge="-0.434" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1443" charge="0.342" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="1444" charge="0.5677" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1445" charge="-0.5881" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1446" charge="0.0713" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1447" charge="0.0985" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1448" charge="-0.0854" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1449" charge="0.0718" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="1450" charge="-0.5232" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1451" charge="0.4422" sigma="0.0" epsilon="0.0"/>
  <Atom type="1452" charge="-0.6318" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="1453" charge="-0.0069" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1454" charge="0.0754" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1455" charge="0.1629" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1456" charge="0.1176" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1457" charge="-0.3691" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1458" charge="0.068" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1459" charge="0.1804" sigma="0.229317330049" epsilon="0.0656888"/>
  <Atom type="1460" charge="-0.0239" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1461" charge="-0.2209" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1462" charge="0.2607" sigma="0.251055258772" epsilon="0.06276"/>
  <Atom type="1463" charge="0.0025" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1464" charge="-0.2269" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1465" charge="0.077" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="1466" charge="0.5194" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1467" charge="-0.5563" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1468" charge="-0.434" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1469" charge="0.342" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="1470" charge="0.5677" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1471" charge="-0.5881" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1472" charge="0.0713" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1473" charge="0.0985" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1474" charge="-0.0854" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1475" charge="0.0718" sigma="0.264953278775" epsilon="0.0656888"/>
  <Atom type="1476" charge="-0.6549" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="1477" charge="0.4396" sigma="0.0" epsilon="0.0"/>
  <Atom type="1478" charge="1.1662" sigma="0.374177461619" epsilon="0.8368"/>
  <Atom type="1479" charge="-0.776" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1480" charge="-0.776" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1481" charge="-0.4989" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1482" charge="0.0558" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1483" charge="0.0679" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1484" charge="0.1065" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1485" charge="0.1174" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1486" charge="-0.3548" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1487" charge="0.0394" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1488" charge="0.2007" sigma="0.229317330049" epsilon="0.0656888"/>
  <Atom type="1489" charge="-0.0251" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1490" charge="0.2006" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1491" charge="0.1553" sigma="0.242146271591" epsilon="0.06276"/>
  <Atom type="1492" charge="-0.6073" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1493" charge="0.0515" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1494" charge="0.7009" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1495" charge="-0.9019" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1496" charge="0.4115" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="1497" charge="-0.7615" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1498" charge="0.5875" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1499" charge="0.0473" sigma="0.242146271591" epsilon="0.06276"/>
  <Atom type="1500" charge="-0.6997" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1501" charge="0.3053" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1502" charge="0.2022" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1503" charge="0.0615" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1504" charge="0.067" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1505" charge="0.0972" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1506" charge="-0.6139" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="1507" charge="0.4186" sigma="0.0" epsilon="0.0"/>
  <Atom type="1508" charge="-0.5246" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1509" charge="1.1662" sigma="0.374177461619" epsilon="0.8368"/>
  <Atom type="1510" charge="-0.776" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1511" charge="-0.776" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1512" charge="-0.4989" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1513" charge="0.0558" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1514" charge="0.0679" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1515" charge="0.1065" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1516" charge="0.1174" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1517" charge="-0.3548" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1518" charge="0.0394" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1519" charge="0.2007" sigma="0.229317330049" epsilon="0.0656888"/>
  <Atom type="1520" charge="-0.0251" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1521" charge="0.2006" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1522" charge="0.1553" sigma="0.242146271591" epsilon="0.06276"/>
  <Atom type="1523" charge="-0.6073" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1524" charge="0.0515" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1525" charge="0.7009" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1526" charge="-0.9019" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1527" charge="0.4115" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="1528" charge="-0.7615" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1529" charge="0.5875" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1530" charge="0.0473" sigma="0.242146271591" epsilon="0.06276"/>
  <Atom type="1531" charge="-0.6997" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1532" charge="0.3053" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1533" charge="0.2022" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1534" charge="0.0615" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1535" charge="0.067" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1536" charge="0.0972" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1537" charge="-0.6139" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="1538" charge="0.4186" sigma="0.0" epsilon="0.0"/>
  <Atom type="1539" charge="-0.6541" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="1540" charge="0.4376" sigma="0.0" epsilon="0.0"/>
  <Atom type="1541" charge="0.4295" sigma="0.0" epsilon="0.0"/>
  <Atom type="1542" charge="-0.6223" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="1543" charge="0.0558" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1544" charge="0.0679" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1545" charge="0.1065" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1546" charge="0.1174" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1547" charge="-0.3548" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1548" charge="0.0394" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1549" charge="0.2007" sigma="0.229317330049" epsilon="0.0656888"/>
  <Atom type="1550" charge="-0.0251" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1551" charge="0.2006" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1552" charge="0.1553" sigma="0.242146271591" epsilon="0.06276"/>
  <Atom type="1553" charge="-0.6073" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1554" charge="0.0515" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1555" charge="0.7009" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1556" charge="-0.9019" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1557" charge="0.4115" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="1558" charge="-0.7615" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1559" charge="0.5875" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1560" charge="0.0473" sigma="0.242146271591" epsilon="0.06276"/>
  <Atom type="1561" charge="-0.6997" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1562" charge="0.3053" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1563" charge="0.2022" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1564" charge="0.0615" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1565" charge="0.067" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1566" charge="0.0972" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1567" charge="-0.6139" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="1568" charge="0.4186" sigma="0.0" epsilon="0.0"/>
  <Atom type="1569" charge="-0.5246" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1570" charge="0.4295" sigma="0.0" epsilon="0.0"/>
  <Atom type="1571" charge="-0.6223" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="1572" charge="0.0558" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1573" charge="0.0679" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1574" charge="0.1065" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1575" charge="0.1174" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1576" charge="-0.3548" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1577" charge="0.0394" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1578" charge="0.2007" sigma="0.229317330049" epsilon="0.0656888"/>
  <Atom type="1579" charge="-0.0251" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1580" charge="0.2006" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1581" charge="0.1553" sigma="0.242146271591" epsilon="0.06276"/>
  <Atom type="1582" charge="-0.6073" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1583" charge="0.0515" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1584" charge="0.7009" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1585" charge="-0.9019" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1586" charge="0.4115" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="1587" charge="-0.7615" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1588" charge="0.5875" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1589" charge="0.0473" sigma="0.242146271591" epsilon="0.06276"/>
  <Atom type="1590" charge="-0.6997" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1591" charge="0.3053" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1592" charge="0.2022" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1593" charge="0.0615" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1594" charge="0.067" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1595" charge="0.0972" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1596" charge="-0.6139" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="1597" charge="0.4186" sigma="0.0" epsilon="0.0"/>
  <Atom type="1598" charge="-0.6541" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="1599" charge="0.4376" sigma="0.0" epsilon="0.0"/>
  <Atom type="1600" charge="1.1662" sigma="0.374177461619" epsilon="0.8368"/>
  <Atom type="1601" charge="-0.776" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1602" charge="-0.776" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1603" charge="-0.4989" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1604" charge="0.0558" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1605" charge="0.0679" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1606" charge="0.1065" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1607" charge="0.1174" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1608" charge="-0.3548" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1609" charge="0.0066" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1610" charge="0.2029" sigma="0.229317330049" epsilon="0.0656888"/>
  <Atom type="1611" charge="-0.0484" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1612" charge="0.0053" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1613" charge="0.1958" sigma="0.251055258772" epsilon="0.06276"/>
  <Atom type="1614" charge="-0.5215" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1615" charge="0.1928" sigma="0.259964245953" epsilon="0.06276"/>
  <Atom type="1616" charge="0.8185" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1617" charge="-0.953" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1618" charge="0.4234" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="1619" charge="-0.7584" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1620" charge="0.7538" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1621" charge="-0.6252" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1622" charge="0.2022" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1623" charge="0.0615" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1624" charge="0.067" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1625" charge="0.0972" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1626" charge="-0.6139" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="1627" charge="0.4186" sigma="0.0" epsilon="0.0"/>
  <Atom type="1628" charge="-0.5246" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1629" charge="1.1662" sigma="0.374177461619" epsilon="0.8368"/>
  <Atom type="1630" charge="-0.776" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1631" charge="-0.776" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1632" charge="-0.4989" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1633" charge="0.0558" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1634" charge="0.0679" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1635" charge="0.1065" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1636" charge="0.1174" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1637" charge="-0.3548" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1638" charge="0.0066" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1639" charge="0.2029" sigma="0.229317330049" epsilon="0.0656888"/>
  <Atom type="1640" charge="-0.0484" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1641" charge="0.0053" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1642" charge="0.1958" sigma="0.251055258772" epsilon="0.06276"/>
  <Atom type="1643" charge="-0.5215" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1644" charge="0.1928" sigma="0.259964245953" epsilon="0.06276"/>
  <Atom type="1645" charge="0.8185" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1646" charge="-0.953" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1647" charge="0.4234" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="1648" charge="-0.7584" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1649" charge="0.7538" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1650" charge="-0.6252" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1651" charge="0.2022" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1652" charge="0.0615" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1653" charge="0.067" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1654" charge="0.0972" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1655" charge="-0.6139" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="1656" charge="0.4186" sigma="0.0" epsilon="0.0"/>
  <Atom type="1657" charge="-0.6541" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="1658" charge="0.4376" sigma="0.0" epsilon="0.0"/>
  <Atom type="1659" charge="0.4295" sigma="0.0" epsilon="0.0"/>
  <Atom type="1660" charge="-0.6223" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="1661" charge="0.0558" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1662" charge="0.0679" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1663" charge="0.1065" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1664" charge="0.1174" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1665" charge="-0.3548" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1666" charge="0.0066" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1667" charge="0.2029" sigma="0.229317330049" epsilon="0.0656888"/>
  <Atom type="1668" charge="-0.0484" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1669" charge="0.0053" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1670" charge="0.1958" sigma="0.251055258772" epsilon="0.06276"/>
  <Atom type="1671" charge="-0.5215" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1672" charge="0.1928" sigma="0.259964245953" epsilon="0.06276"/>
  <Atom type="1673" charge="0.8185" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1674" charge="-0.953" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1675" charge="0.4234" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="1676" charge="-0.7584" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1677" charge="0.7538" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1678" charge="-0.6252" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1679" charge="0.2022" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1680" charge="0.0615" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1681" charge="0.067" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1682" charge="0.0972" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1683" charge="-0.6139" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="1684" charge="0.4186" sigma="0.0" epsilon="0.0"/>
  <Atom type="1685" charge="-0.5246" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1686" charge="0.4295" sigma="0.0" epsilon="0.0"/>
  <Atom type="1687" charge="-0.6223" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="1688" charge="0.0558" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1689" charge="0.0679" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1690" charge="0.1065" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1691" charge="0.1174" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1692" charge="-0.3548" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1693" charge="0.0066" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1694" charge="0.2029" sigma="0.229317330049" epsilon="0.0656888"/>
  <Atom type="1695" charge="-0.0484" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1696" charge="0.0053" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1697" charge="0.1958" sigma="0.251055258772" epsilon="0.06276"/>
  <Atom type="1698" charge="-0.5215" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1699" charge="0.1928" sigma="0.259964245953" epsilon="0.06276"/>
  <Atom type="1700" charge="0.8185" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1701" charge="-0.953" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1702" charge="0.4234" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="1703" charge="-0.7584" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1704" charge="0.7538" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1705" charge="-0.6252" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1706" charge="0.2022" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1707" charge="0.0615" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1708" charge="0.067" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1709" charge="0.0972" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1710" charge="-0.6139" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="1711" charge="0.4186" sigma="0.0" epsilon="0.0"/>
  <Atom type="1712" charge="-0.6541" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="1713" charge="0.4376" sigma="0.0" epsilon="0.0"/>
  <Atom type="1714" charge="1.1662" sigma="0.374177461619" epsilon="0.8368"/>
  <Atom type="1715" charge="-0.776" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1716" charge="-0.776" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1717" charge="-0.4989" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1718" charge="0.0558" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1719" charge="0.0679" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1720" charge="0.1065" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1721" charge="0.1174" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1722" charge="-0.3548" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1723" charge="0.0191" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1724" charge="0.2006" sigma="0.229317330049" epsilon="0.0656888"/>
  <Atom type="1725" charge="0.0492" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1726" charge="0.1374" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1727" charge="0.164" sigma="0.242146271591" epsilon="0.06276"/>
  <Atom type="1728" charge="-0.5709" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1729" charge="0.1744" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1730" charge="0.477" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1731" charge="-0.5597" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1732" charge="-0.4787" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1733" charge="0.3424" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="1734" charge="0.7657" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1735" charge="-0.9672" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1736" charge="0.4364" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="1737" charge="-0.6323" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1738" charge="0.1222" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1739" charge="0.2022" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1740" charge="0.0615" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1741" charge="0.067" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1742" charge="0.0972" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1743" charge="-0.6139" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="1744" charge="0.4186" sigma="0.0" epsilon="0.0"/>
  <Atom type="1745" charge="-0.5246" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1746" charge="1.1662" sigma="0.374177461619" epsilon="0.8368"/>
  <Atom type="1747" charge="-0.776" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1748" charge="-0.776" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1749" charge="-0.4989" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1750" charge="0.0558" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1751" charge="0.0679" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1752" charge="0.1065" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1753" charge="0.1174" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1754" charge="-0.3548" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1755" charge="0.0191" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1756" charge="0.2006" sigma="0.229317330049" epsilon="0.0656888"/>
  <Atom type="1757" charge="0.0492" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1758" charge="0.1374" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1759" charge="0.164" sigma="0.242146271591" epsilon="0.06276"/>
  <Atom type="1760" charge="-0.5709" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1761" charge="0.1744" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1762" charge="0.477" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1763" charge="-0.5597" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1764" charge="-0.4787" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1765" charge="0.3424" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="1766" charge="0.7657" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1767" charge="-0.9672" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1768" charge="0.4364" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="1769" charge="-0.6323" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1770" charge="0.1222" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1771" charge="0.2022" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1772" charge="0.0615" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1773" charge="0.067" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1774" charge="0.0972" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1775" charge="-0.6139" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="1776" charge="0.4186" sigma="0.0" epsilon="0.0"/>
  <Atom type="1777" charge="-0.6541" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="1778" charge="0.4376" sigma="0.0" epsilon="0.0"/>
  <Atom type="1779" charge="0.4295" sigma="0.0" epsilon="0.0"/>
  <Atom type="1780" charge="-0.6223" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="1781" charge="0.0558" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1782" charge="0.0679" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1783" charge="0.1065" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1784" charge="0.1174" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1785" charge="-0.3548" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1786" charge="0.0191" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1787" charge="0.2006" sigma="0.229317330049" epsilon="0.0656888"/>
  <Atom type="1788" charge="0.0492" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1789" charge="0.1374" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1790" charge="0.164" sigma="0.242146271591" epsilon="0.06276"/>
  <Atom type="1791" charge="-0.5709" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1792" charge="0.1744" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1793" charge="0.477" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1794" charge="-0.5597" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1795" charge="-0.4787" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1796" charge="0.3424" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="1797" charge="0.7657" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1798" charge="-0.9672" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1799" charge="0.4364" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="1800" charge="-0.6323" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1801" charge="0.1222" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1802" charge="0.2022" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1803" charge="0.0615" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1804" charge="0.067" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1805" charge="0.0972" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1806" charge="-0.6139" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="1807" charge="0.4186" sigma="0.0" epsilon="0.0"/>
  <Atom type="1808" charge="-0.5246" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1809" charge="0.4295" sigma="0.0" epsilon="0.0"/>
  <Atom type="1810" charge="-0.6223" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="1811" charge="0.0558" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1812" charge="0.0679" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1813" charge="0.1065" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1814" charge="0.1174" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1815" charge="-0.3548" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1816" charge="0.0191" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1817" charge="0.2006" sigma="0.229317330049" epsilon="0.0656888"/>
  <Atom type="1818" charge="0.0492" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1819" charge="0.1374" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1820" charge="0.164" sigma="0.242146271591" epsilon="0.06276"/>
  <Atom type="1821" charge="-0.5709" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1822" charge="0.1744" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1823" charge="0.477" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1824" charge="-0.5597" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1825" charge="-0.4787" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1826" charge="0.3424" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="1827" charge="0.7657" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1828" charge="-0.9672" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1829" charge="0.4364" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="1830" charge="-0.6323" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1831" charge="0.1222" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1832" charge="0.2022" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1833" charge="0.0615" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1834" charge="0.067" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1835" charge="0.0972" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1836" charge="-0.6139" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="1837" charge="0.4186" sigma="0.0" epsilon="0.0"/>
  <Atom type="1838" charge="-0.6541" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="1839" charge="0.4376" sigma="0.0" epsilon="0.0"/>
  <Atom type="1840" charge="1.1662" sigma="0.374177461619" epsilon="0.8368"/>
  <Atom type="1841" charge="-0.776" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1842" charge="-0.776" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1843" charge="-0.4989" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1844" charge="0.0558" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1845" charge="0.0679" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1846" charge="0.1065" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1847" charge="0.1174" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1848" charge="-0.3548" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1849" charge="0.0674" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1850" charge="0.1824" sigma="0.229317330049" epsilon="0.0656888"/>
  <Atom type="1851" charge="0.0418" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1852" charge="-0.1126" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1853" charge="0.2188" sigma="0.251055258772" epsilon="0.06276"/>
  <Atom type="1854" charge="-0.3635" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1855" charge="0.1811" sigma="0.259964245953" epsilon="0.06276"/>
  <Atom type="1856" charge="0.5952" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1857" charge="-0.5761" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1858" charge="-0.3549" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1859" charge="0.3154" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="1860" charge="0.4687" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1861" charge="-0.5477" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1862" charge="0.2022" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1863" charge="0.0615" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1864" charge="0.067" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1865" charge="0.0972" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1866" charge="-0.6139" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="1867" charge="0.4186" sigma="0.0" epsilon="0.0"/>
  <Atom type="1868" charge="-0.5246" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1869" charge="1.1662" sigma="0.374177461619" epsilon="0.8368"/>
  <Atom type="1870" charge="-0.776" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1871" charge="-0.776" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1872" charge="-0.4989" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1873" charge="0.0558" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1874" charge="0.0679" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1875" charge="0.1065" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1876" charge="0.1174" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1877" charge="-0.3548" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1878" charge="0.0674" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1879" charge="0.1824" sigma="0.229317330049" epsilon="0.0656888"/>
  <Atom type="1880" charge="0.0418" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1881" charge="-0.1126" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1882" charge="0.2188" sigma="0.251055258772" epsilon="0.06276"/>
  <Atom type="1883" charge="-0.3635" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1884" charge="0.1811" sigma="0.259964245953" epsilon="0.06276"/>
  <Atom type="1885" charge="0.5952" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1886" charge="-0.5761" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1887" charge="-0.3549" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1888" charge="0.3154" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="1889" charge="0.4687" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1890" charge="-0.5477" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1891" charge="0.2022" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1892" charge="0.0615" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1893" charge="0.067" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1894" charge="0.0972" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1895" charge="-0.6139" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="1896" charge="0.4186" sigma="0.0" epsilon="0.0"/>
  <Atom type="1897" charge="-0.6541" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="1898" charge="0.4376" sigma="0.0" epsilon="0.0"/>
  <Atom type="1899" charge="0.4295" sigma="0.0" epsilon="0.0"/>
  <Atom type="1900" charge="-0.6223" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="1901" charge="0.0558" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1902" charge="0.0679" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1903" charge="0.1065" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1904" charge="0.1174" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1905" charge="-0.3548" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1906" charge="0.0674" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1907" charge="0.1824" sigma="0.229317330049" epsilon="0.0656888"/>
  <Atom type="1908" charge="0.0418" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1909" charge="-0.1126" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1910" charge="0.2188" sigma="0.251055258772" epsilon="0.06276"/>
  <Atom type="1911" charge="-0.3635" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1912" charge="0.1811" sigma="0.259964245953" epsilon="0.06276"/>
  <Atom type="1913" charge="0.5952" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1914" charge="-0.5761" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1915" charge="-0.3549" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1916" charge="0.3154" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="1917" charge="0.4687" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1918" charge="-0.5477" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1919" charge="0.2022" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1920" charge="0.0615" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1921" charge="0.067" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1922" charge="0.0972" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1923" charge="-0.6139" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="1924" charge="0.4186" sigma="0.0" epsilon="0.0"/>
  <Atom type="1925" charge="-0.5246" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1926" charge="0.4295" sigma="0.0" epsilon="0.0"/>
  <Atom type="1927" charge="-0.6223" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="1928" charge="0.0558" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1929" charge="0.0679" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1930" charge="0.1065" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1931" charge="0.1174" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1932" charge="-0.3548" sigma="0.300001234347" epsilon="0.71128"/>
  <Atom type="1933" charge="0.0674" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1934" charge="0.1824" sigma="0.229317330049" epsilon="0.0656888"/>
  <Atom type="1935" charge="0.0418" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1936" charge="-0.1126" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1937" charge="0.2188" sigma="0.251055258772" epsilon="0.06276"/>
  <Atom type="1938" charge="-0.3635" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1939" charge="0.1811" sigma="0.259964245953" epsilon="0.06276"/>
  <Atom type="1940" charge="0.5952" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1941" charge="-0.5761" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1942" charge="-0.3549" sigma="0.324999852378" epsilon="0.71128"/>
  <Atom type="1943" charge="0.3154" sigma="0.106907846177" epsilon="0.0656888"/>
  <Atom type="1944" charge="0.4687" sigma="0.339966950842" epsilon="0.359824"/>
  <Atom type="1945" charge="-0.5477" sigma="0.295992190115" epsilon="0.87864"/>
  <Atom type="1946" charge="0.2022" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1947" charge="0.0615" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1948" charge="0.067" sigma="0.339966950842" epsilon="0.4577296"/>
  <Atom type="1949" charge="0.0972" sigma="0.247135304412" epsilon="0.0656888"/>
  <Atom type="1950" charge="-0.6139" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="1951" charge="0.4186" sigma="0.0" epsilon="0.0"/>
  <Atom type="1952" charge="-0.6541" sigma="0.306647338784" epsilon="0.8803136"/>
  <Atom type="1953" charge="0.4376" sigma="0.0" epsilon="0.0"/>
  <Atom type="1954" charge="-1.0" sigma="0.440103966761" epsilon="0.4184"/>
  <Atom type="1955" charge="1.0" sigma="0.604920229617" epsilon="0.0003372304"/>
  <Atom type="1956" charge="1.0" sigma="0.473601758563" epsilon="0.001372352"/>
  <Atom type="1957" charge="1.0" sigma="0.202590368505" epsilon="0.0765672"/>
  <Atom type="1958" charge="2.0" sigma="0.1412252648" epsilon="3.7434248"/>
  <Atom type="1959" charge="1.0" sigma="0.332839761097" epsilon="0.01158968"/>
  <Atom type="1960" charge="1.0" sigma="0.526699322165" epsilon="0.00071128"/>
 </NonbondedForce>''')
        print('XML MADE!')


protein = QUBEMAKER('ionized.pdb', PAR_file='QUBE_FF_openmm.par')