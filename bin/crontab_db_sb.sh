

# every six months updates the databases
0 0 1 */6 * "d:/projects/databases/bin/create_fasta-sb.sh" &> "d:/projects/databases/logs/create_fasta-sb.log"
#0 0 4 */1 * "d:/projects/databases/bin/create_fasta-sb.sh" &> "d:/projects/databases/logs/create_fasta-sb.log"

