import lomap

# If multiple toolkits are present RDK will be the default one
db_rdk = lomap.DBMolecules("./tests/data")

# Checking RDK mol
print(db_rdk[0].molecule)

# Switching Toolkit

try:
    lomap.toolkits.set_default("OE")

    db_oe = lomap.DBMolecules("./tests/data")

    # Checking OE mol
    print(db_oe[0].molecule)
except:
    pass
