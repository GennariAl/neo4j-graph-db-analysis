### INDICE ###
# RICERCA WCC
# 2. runno il Weakly Connected Component algorithm per identificare le community (componenti connesse)
# 3. creo la lista wcc, che contiene un dizionario per ogni community (id = , subgs = , samples = )

import logging
import statistics
from neo4j import GraphDatabase
from neo4j.exceptions import Neo4jError
from difflib import SequenceMatcher
from Bio import Align

uri = "bolt://localhost:7687"
user = "neo4j"
password = "Cambiami2020!"

def get_driver():
    return GraphDatabase.driver
get_driver()
driver = GraphDatabase.driver(uri, auth=(user, password))
driver.verify_connectivity()

# 1. projection creation
def projection(tx):
	result = tx.run("""call gds.graph.project('subg-interactions', 'subg', 'gtris');""")
	return [ record for record in result ]

# with driver.session() as session:
#     result = session.write_transaction(projection)
#     for record in result:
#       print(record)
#     session.close()

# 2. weakly connected components indentification
def get_community(tx):
	result = tx.run("""call gds.wcc.stream('subg-interactions', {})
					yield componentId, nodeId 
					return gds.util.asNode(nodeId).clu_id, componentId 
					order by componentId;""")
	return [ record for record in result ]

subg_list = []
community_list = []		# two equal long lists, overlapping them I get the clu_id and componentId
with driver.session() as session:
	result = session.read_transaction(get_community)
	for record in result:
		subg_list.append(record["gds.util.asNode(nodeId).clu_id"])
		community_list.append(record["componentId"])
	session.close()

# 3. I create the WCC list, which contains a dictionary for each community (id = , subgs = , samples = )
# to count the subgs in each community I just count the occurrence of each componentId in the community_list

community_dict = {}		# dictionary that collects communities and the number of subg's that have
for i in range(0, len(community_list)):
	if community_list[i] not in community_dict.keys():
		community_dict[community_list[i]] = 1
	else:
		community_dict[community_list[i]] = community_dict[community_list[i]] + 1

community_dict = dict(sorted(community_dict.items(), key = lambda item: item[1], reverse=True))






# funzione per estrarre il sample di un subg
def samples_extractor(tx, s):
	result = tx.run("""match (s:sample)-[r:sample2subg]-(u:subg)
					where u.clu_id = $subg
					return s.CompleteAmplificationID;""", subg = s)
	return [ record for record in result ]

wcc = []   	# list that will contain the dictionaries of each community with id, number of nodes, average centrality 
count = 0
for c in community_dict.keys():

	if community_dict[c] > 600:		# FILTER FOR CONSIDERING A PORTION OF THE TOTAL COMPONENTS

		print(f'\n> I CONSIDER THE COMMUNITY WITH CON {community_dict[c]} NODES')
		samples_in_comm = []
		subg_in_comm = []
		comm_dict_completo = {}
		for i in range(0, len(community_list)):		# I iterate over the community_list to isolate the one with Id = c
			if community_list[i] == c:
				subg_in_comm.append(subg_list[i])

				with driver.session() as session:
					result = session.write_transaction(samples_extractor, subg_list[i])
					for record in result:
						if record["s.CompleteAmplificationID"] not in samples_in_comm:
							samples_in_comm.append(record["s.CompleteAmplificationID"])
				session.close()

		print(f'> numero di subg raccolti: {len(subg_in_comm)}')
		print(f'> numero di sample raccolti: {len(samples_in_comm)}')
		comm_dict_completo['id'] = c
		comm_dict_completo['subgs'] = len(subg_in_comm)
		comm_dict_completo['samples'] = len(samples_in_comm)
		wcc.append(comm_dict_completo)

		count += 1
		print(count)

# function to compute similarity score
def mysimilarity(a, b):
	a = a[:50]				# considero solo le basi iniziali, perchè sono quelle veramente responsabili dell'allineamento sul genoma
	b = b[:50]
	match_scores = []		# contiene gli scores di similarity per ogni caso

	# caso -2
	new_a = a[2:]		# seq a shifta a sinistra di due basi, quindi inizia da 2
	#print(f'{new_a}\n{b}\n')
	lennew_a = len(new_a)
	lenb = len(b)
	match = 0
	for i in range(0, min(lennew_a, lenb)):	
		if new_a[i] == b[i]:
			match += 1
	match_scores.append(match)

	# caso -1
	new_a = a[1:]		# seq a shifta a sinistra di una base, quindi inizia da 1
	#print(f'{new_a}\n{b}\n')
	lennew_a = len(new_a)
	lenb = len(b)
	match = 0
	for i in range(0, min(lennew_a, lenb)):	
		if new_a[i] == b[i]:
			match += 1
	match_scores.append(match)

	# caso 0
	#print(f'{a}\n{b}\n')
	lena = len(a)
	lenb = len(b)
	match = 0
	for i in range(0, min(lena, lenb)):	
		if a[i] == b[i]:
			match += 1
	match_scores.append(match)
	
	# caso +1
	new_b = b[1:]		# seq b shifta a sinistra di una base, quindi inizia da 1
	#print(f'{a}\n{new_b}\n')
	lena = len(a)
	lennewb = len(new_b)
	match = 0
	for i in range(0, min(lena, lennewb)):	
		if a[i] == new_b[i]:
			match += 1
	match_scores.append(match)
	    
	# caso +2
	new_b = b[2:]		# seq b shifta a sinistra di due basi, quindi inizia da 2
	#print(f'{a}\n{new_b}')
	lena = len(a)
	lennewb = len(new_b)
	match = 0
	for i in range(0, min(lena, lennewb)):	
		if a[i] == new_b[i]:
			match += 1
	match_scores.append(match)

	#print(match_scores)
	#print(f'{max(match_scores)} / {min([lena,lenb])}')
	return max(match_scores)/min([lena,lenb])

def cons_seq_extractor(tx, clu_id):
	result = tx.run("""match (u:subg)
					where u.clu_id = $clu_id
					return u.cons_seq;""", clu_id = clu_id)
	return [ record for record in result ]

def connection_extractor(tx, clu_id):
	result = tx.run("""match (u:subg)-[r:gtris]-(b:subg)
					where u.clu_id = $clu_id
					return r.case, b.clu_id, b.cons_seq;""", clu_id = clu_id)
	return [ record for record in result ]

worse_cases = []
for c in wcc:
	subg_in_community = []
	simi_scores_in_component = []
	# I create the list of all nodes subg in the component c
	for i in range(0, len(community_list)):
		if community_list[i] == c['id']:
			subg_in_community.append(subg_list[i])

	# now for each subg node I extract from the graphDB the consensus sequence
	dizionario_subg_cons_seq = {}
	for i in subg_in_community:
		with driver.session() as session:
			result = session.read_transaction(cons_seq_extractor, i)
		for record in result:
			dizionario_subg_cons_seq[i] = record["u.cons_seq"]

	# For each subg node I calculate the similarity with all other subg that it connects
	print(f'> inizio ad calcolare le similarity nel component con {c["subgs"]} subg')
	for k in dizionario_subg_cons_seq.keys():	
		scores_della_subg = []	
		with driver.session() as session:
			result = session.read_transaction(connection_extractor, k)		# I get all outgoing connections from k
			for record in result:											# for each connection I calculate the similarity score
				similarità = mysimilarity(dizionario_subg_cons_seq[k], record['b.cons_seq'])
				scores_della_subg.append(similarità)
		simi_scores_in_component.append(statistics.mean(scores_della_subg))

	# file .csv collecting for each component, the average score of all subg nodes
	filename = 'SIMILARITY_'+str(c["subgs"])+'.csv'
	with open(filename, 'w') as f:
		f.write('similarità\n')
		for i in simi_scores_in_component:
			f.write(f'{i}\n')
