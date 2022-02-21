import pandas as pd
from sqlalchemy import create_engine, MetaData, Table, select
from sqlalchemy.engine.url import URL


config = {"drivername": "postgresql",
          "host": "localhost",
          "port": "5432",
          "username": "chen1i6c04",
          "password": "tp61i6c04"}


def connect_database(database):
    config['database'] = database
    engine = create_engine(URL(**config))
    con = engine.connect()
    metadata = MetaData()
    return engine, con, metadata


def load_table(database, table_name):
    engine, con, metadata = connect_database(database)
    query = select('*', Table(table_name, metadata))
    result = pd.read_sql(query, con)
    con.close()
    engine.dispose()
    return result
