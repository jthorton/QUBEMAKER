import xml.etree.ElementTree as ET

# Now parse the xml file for the rest of the data
inputXML_file = 'test.xml'
inXML = ET.parse(inputXML_file)
in_root = inXML.getroot()

atom_names2 = ['CAY', 'HY1', 'HY2', 'HY3', 'CY', 'OY', 'N', 'HN', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'HG', 'CD1', 'HD11', 'HD12', 'HD13', 'CD2', 'HD21', 'HD22', 'HD23', 'C', 'O', 'NT', 'HNT', 'CAT', 'HT1', 'HT2', 'HT3']


atom_map = {'NT': 'N706', 'HNT': 'H707', 'CAT': 'C708', 'HT1': 'H709', 'HT2': 'H709', 'HT3': 'H709', 'CAY':'C711', 'HY1': 'H710', 'HY2': 'H710', 'HY3': 'H710', 'CY': 'C712', 'OY': 'O713', 'N': 'N200', 'HN': 'H201', 'CA': 'C202', 'HA' : 'H203', 'CB': 'C204', 'HB1': 'H205', 'HB2': 'H205', 'CG': 'C206', 'HG': 'H207', 'CD1': 'C208', 'HD11': 'H209', 'HD12': 'H209', 'HD13': 'H209', 'CD2': 'C210', 'HD21': 'H211', 'HD22': 'H211', 'HD23': 'H211', 'C': 'C212', 'O': 'O213'}


exception_dict1 = {}
for Exception in in_root.iter('Exception'):
    if Exception.get("q") != "0":
        # print(f'{atom_names2[int(Exception.get("p1"))]}-{atom_names2[int(Exception.get("p2"))]} = {Exception.get("q")}')
        exception_dict1[(atom_names2[int(Exception.get("p1"))], atom_names2[int(Exception.get("p2"))])] = Exception.get("q")




# angles1 = 0
# angle_list_1 =[]
#
# for Angle in in_root.iter('Angle'):
#     try :
#         angles1 +=1
#         angle_list_1.append((atom_names2[int(Angle.get("p1"))], atom_names2[int(Angle.get("p2"))], atom_names2[int(Angle.get("p3"))]))
#         print(f'<Angle class1="{atom_names2[int(Angle.get("p1"))]}" class2="{atom_names2[int(Angle.get("p2"))]}" class3="{atom_names2[int(Angle.get("p3"))]}" angle="{Angle.get("a")}" k="{Angle.get("k")}"/>')
#     except KeyError:
#          print('missing')

# torsions = 0
# torsion_list_1 = []
# for torsion in in_root.iter('Torsion'):
#     try:
#         torsions += 1
#         torsion_list_1.append((atom_names2[int(torsion.get("p1"))], atom_names2[int(torsion.get("p2"))], atom_names2[int(torsion.get("p3"))], atom_names2[int(torsion.get("p4"))]))
#         print(f'<Torsion class1="{atom_map[atom_names2[int(torsion.get("p1"))]]}" class2="{atom_map[atom_names2[int(torsion.get("p2"))]]}" class3="{atom_map[atom_names2[int(torsion.get("p3"))]]}" class4="{atom_map[atom_names2[int(torsion.get("p4"))]]}" periodicity="{torsion.get("periodicity")}" phase="{torsion.get("phase")}" k="{torsion.get("k")}"/>')
#     except KeyError:
#         pass
#
# print(torsions)

#angle_list_2 =[]
#
inputXML_file = 'serialized_system.xml'
# inputXML_file = 'test.xml'
inXML = ET.parse(inputXML_file)
in_root = inXML.getroot()

atom_names = ['NT', 'HNT', 'CAT', 'HT1', 'HT2', 'HT3', 'CAY', 'HY1', 'HY2', 'HY3', 'CY', 'OY', 'N', 'HN', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'HG', 'CD1', 'HD11', 'HD12', 'HD13', 'CD2', 'HD21', 'HD22', 'HD23', 'C', 'O']




exception_dict2 = {}
for Exception in in_root.iter('Exception'):
    if Exception.get("q") != "0":
        # print(f'{atom_names[int(Exception.get("p1"))]}-{atom_names[int(Exception.get("p2"))]} = {Exception.get("q")}')
        exception_dict2[(atom_names[int(Exception.get("p1"))], atom_names[int(Exception.get("p2"))])] = Exception.get("q")


# now loop over all of the exceptions
print(exception_dict1)
print(exception_dict2)
for key in exception_dict2.keys():
    try:
        if exception_dict2[key] != exception_dict1[key]:
            print(key, exception_dict1[key], exception_dict2[key])
    except KeyError:
        if exception_dict2[key] != exception_dict1[key[::-1]]:
            print(key, exception_dict1[key[::-1]], exception_dict2[key])



#
# angles2 = 0
# for Angle in in_root.iter('Angle'):
#     try :
#         angles2 +=1
#         angle_list_2.append((atom_names[int(Angle.get("p1"))], atom_names[int(Angle.get("p2"))], atom_names[int(Angle.get("p3"))]))
#         print(f'<Angle class1="{atom_names[int(Angle.get("p1"))]}" class2="{atom_names[int(Angle.get("p2"))]}" class3="{atom_names[int(Angle.get("p3"))]}" angle="{Angle.get("a")}" k="{Angle.get("k")}"/>')
#     except KeyError:
#         print('missing')
#
#
# print(f'angles test = {angles1}')
# print(f'angles serialized = {angles2}')
#
# print(f'test angles = {angle_list_1}')
# print(f'serialized angles = {angle_list_2}')
# # find the missing angle
# for angle in angle_list_2:
#     if angle not in angle_list_1:
#         if angle[::-1] not in angle_list_1:
#             print(atom_map[angle[0]], atom_map[angle[1]], atom_map[angle[2]])

# torsions = 0
# torsion_list_2 = []
# for torsion in in_root.iter('Torsion'):
#     try:
#         torsions += 1
#         torsion_list_2.append((atom_names[int(torsion.get("p1"))], atom_names[int(torsion.get("p2"))],
#                                atom_names[int(torsion.get("p3"))], atom_names[int(torsion.get("p4"))]))
#         print(f'<Torsion class1="{atom_map[atom_names[int(torsion.get("p1"))]]}" class2="{atom_map[atom_names[int(torsion.get("p2"))]]}" class3="{atom_map[atom_names[int(torsion.get("p3"))]]}" class4="{atom_map[atom_names[int(torsion.get("p4"))]]}" periodicity="{torsion.get("periodicity")}" phase="{torsion.get("phase")}" k="{torsion.get("k")}"/>')
#     except KeyError:
#         pass
#
# print(torsions)
#
# # find the missing torsions
# for torsion in torsion_list_2:
#     if torsion not in torsion_list_1:
#         if torsion[::-1] not in torsion_list_1:
#             print(atom_map[torsion[0]], atom_map[torsion[1]], atom_map[torsion[2]], atom_map[torsion[3]])