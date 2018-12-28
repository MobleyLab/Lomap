import lomap

# If multiple toolkits are present OE will be the default one
db_oe = lomap.DBMolecules("./tests/data")

# Checking OE mol
print(db_oe[0].molecule)

# Switching Toolkit
lomap.toolkits.set_default("RDK")

db_rdk = lomap.DBMolecules("./tests/data")

# Checking RDK mol
print(db_rdk[0].molecule)
