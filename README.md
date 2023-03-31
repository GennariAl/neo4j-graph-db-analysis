# neo4j-graph-db-analysis-in-Python
Each analysis on the graph DB has been carried out through Python scripts which queries the DB and processes the results

## 1. Weakly Connected Components identification
To connect to and query Neo4j from within a Python application, it is used the Neo4j Python Driver: a thread-safe, application-wide fixture from which all Neo4j interaction derives. In order to query the database, through the Driver it is required to open Sessions. A
session is a container for a sequence of transactions. Sessions borrow connections from a pool as required and are considered lightweight and disposable. Through a Session, it is possible to run one or more Transactions, each transaction comprises a unit of work performed against a database. It is treated in a coherent and reliable way, independent of other transactions. In Neo4j Graph Data Science Library (GDS), algorithms can be executed on a projection (in this analysis called "subg-interactions"): a named graph that
has been filtered based on its node labels and relationship types, that filtered graph only exists during the execution of the algorithm.
The complete description of the algorithm is available here:


## 2. Sequence similarity score computation
To assess how much two sequences are identical, a similarity score is computed as the number of matches over the maximum number of matches possible, which correspond to the length of the shortest sequence, and it goes from a minimum of 0 to a maximum of 1. In contrast to the alignment, the similarity score com- putation starts from the first position of the two sequences and no gaps are added. 
The function I implemented takes each position of the two sequences and if, the two corresponding nucleotides are the same, increase a ’match’ variable by 1. Then moves to the next position and repeat the process, until it reach the shortest of the two sequences. The ratio of number of matches, corresponding to the variable ’match’, over the maximum number of matches (corresponding to the shortest sequence length) is then returned as similarity score. The complete description of the algorithm, and its use on the WCC algorithm results is available here:


