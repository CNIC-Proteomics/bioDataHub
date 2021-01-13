
# test
* * * * * env > "d:/tmp/env.out"
* * * * * echo "test 1" >> "d:/tmp/test.log"

# every six months updates the databases
0 0 1 */6 * "d:/projects/databases/bin/create_db_sb.sh" &> "d:/projects/databases/logs/create_db_sb.log"
#0 0 4 */1 * "d:/projects/databases/bin/create_db_sb.sh" &> "d:/projects/databases/logs/create_db_sb.log"

