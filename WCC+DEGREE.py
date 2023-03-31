import logging
import statistics
from neo4j import GraphDatabase
from neo4j.exceptions import Neo4jError

uri = "bolt://localhost:7687"
user = "neo4j"
password = "Cambiami2020!"

# driver
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
def get_comm(tx):
	result = tx.run("""call gds.wcc.stream('subg-interactions', {})
					yield componentId, nodeId 
					return gds.util.asNode(nodeId).clu_id, componentId 
					order by componentId;""")
	return [ record for record in result ]

subg_list = []
community_list = []		# two equal long lists, overlapping them I get the clu_id and componentId
with driver.session() as session:
	result = session.read_transaction(get_comm)
	for record in result:
		subg_list.append(record["gds.util.asNode(nodeId).clu_id"])
		community_list.append(record["componentId"])
	session.close()

print(f'> len comm_list: {len(community_list)}')
print(f'> len subg_list: {len(subg_list)}')

# 3. I create the WCC list, which contains a dictionary for each community (id = , subgs = , samples = )
# to count the subgs in each community I just count the occurrence of each componentId in the community_list

community_dict = {}		# dictionary that collects communities and the number of subg's that have
for i in range(0, len(community_list)):
	if community_list[i] not in community_dict.keys():
		community_dict[community_list[i]] = 1
	else:
		community_dict[community_list[i]] = community_dict[community_list[i]] + 1

community_dict = dict(sorted(community_dict.items(), key = lambda item: item[1], reverse=True))

# 4. I calculate the centrality of all nodes subg
print(f'> calculate the centralities')
def centrality(tx):
	result = tx.run(	"""call gds.degree.stream('subg-interactions', {orientation:"UNDIRECTED"}) yield nodeId, score return gds.util.asNode(nodeId).clu_id, score;""")
	return [ record for record in result ]

subg_centrality = {}
with driver.session() as session:
	result = session.read_transaction(centrality)
	for record in result:
		subg_centrality[record["gds.util.asNode(nodeId).clu_id"]] = record["score"]
	session.close()

wcc_degree = []   	# list that will contain the dictionaries of each community with id, number of nodes, average centrality 
count = 0

for c in community_dict.keys():
	if community_dict[c] > 1:		# FILTER FOR CONSIDERING A PORTION OF THE TOTAL COMPONENTS

		print(f'> I CONSIDER THE COMMUNITY WITH CON {community_dict[c]} NODES')
		samples_in_comm = []
		subg_in_community = []

		for i in range(0, len(community_list)):		# I iterate over the community_list to isolate the one with Id = c
			if community_list[i] == c:
				subg_in_community.append(subg_list[i])

		# I extract only the centrality of subg nodes in the community
		chosen_centrality = []						# list of all centralities of all subg's in the component
		for i in subg_in_community:					# loop over the list of subg nodes in the community
			if i in subg_centrality.keys():			# if the subg has a centrality score
				chosen_centrality.append(subg_centrality[i])
			else:
				print('\tSUBG HAS NO CENTRALITY SCORE')

		comm_dict_completo = {}
		comm_dict_completo['id'] = c
		comm_dict_completo['subgs'] = len(subg_in_community)
		comm_dict_completo['centrality'] = statistics.mean(chosen_centrality)
		wcc_degree.append(comm_dict_completo)

		count += 1
		print(count)

# write a .csv file with all the components and their values (contained in the dictionary)
filename = str('WCC_DEGREE.csv')
with open(filename, 'w') as f:
	f.write(f'WCC,SUBG NODES,AVG DEGREE\n')
	for i in wcc_degree:
			f.write(f'{i["id"]},{i["subgs"]},{i["centrality"]}\n')
			