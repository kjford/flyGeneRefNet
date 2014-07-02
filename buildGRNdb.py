'''
buildGRNdb.py

Construct GeneRefNet MySQL database from Flybase data (flybase.org)

Copyright Kevin Ford (2014)

dependencies:
Packages you'll likely need to find on your own:
psycopg2 (http://initd.org/psycopg/)
MySQLdb (http://mysql-python.sourceforge.net/MySQLdb.html)
Packages typically installed with a Canopy/Enthought install of Python:
pandas
numpy

Note on usage:
For running your own MySQL database you will need full privileges to a database
(i.e. create, insert, select, etc.)
The code below uses a local database called 'flydb' with authentication given below.
You will need to create the database and user in MySQL first.
This is obviously not very secure since the authentication here is in plain text.
Once built, you can revoke all privileges except SELECT.

Tables:
pub: Pub id, title, reference, abstract
flygenes: Unique Symb ID, FBgnID, symbol
generefs: unique generefID, Pub id, FBgnID

'''
import psycopg2 as psdb
import MySQLdb as mdb
import pandas as pd
import numpy as np

# Flybase authentication:
authent_flybase = {'host':'flybase.org',
           'user':'flybase',
           'dbname':'flybase'}

# local authentication
# change to reflect your login, password and schema name:
authent_local = {'host':'localhost',
           'user':'pyuser',
           'passwd':'testpass',
           'db':'flydb'}

'''
Functions for pulling data from Flybase postres database
for schema organization see: http://gmod.org/wiki/FlyBase_Field_Mapping_Tables
'''

# pull unique ids and gene symbols from Flybase
def pullFBgns(auth=authent_flybase):
    qry='''
    SELECT feature.uniquename, feature.name
    FROM feature, cvterm, organism o
    WHERE cvterm.name='gene' AND
    feature.type_id=cvterm.cvterm_id AND
    feature.is_obsolete='f' AND
    feature.uniquename LIKE 'FBgn%' AND
    feature.organism_id=o.organism_id AND
    o.species='melanogaster'
    '''
    con = psdb.connect(**auth)
    with con:
        print('Getting IDs and gene symbols from Flybase')
        print('Connecting to Flybase')
        df=pd.io.sql.read_frame(qry,con)
        print('Found %i records'%df.shape[0])
    return df

# pull all references that are papers from flybase with titles and abstracts
def pullFBrefs(auth=authent_flybase):
    qry='''
    SELECT A.pub_id, A.title, A.miniref, B.abst
    FROM (
        SELECT pub.pub_id pub_id, pub.title title, pub.miniref miniref
        FROM pub, cvterm cv1
        WHERE pub.type_id=cv1.cvterm_id AND
        cv1.name='paper') A
    LEFT OUTER JOIN (
        SELECT pp.VALUE abst, pp.pub_id pid
        FROM cvterm cv2, pubprop pp
        WHERE cv2.name='pubmed_abstract' AND
        pp.type_id=cv2.cvterm_id) B
    ON (A.pub_id=B.pid)
    '''
    con = psdb.connect(**auth)
    with con:
        print('Getting references with abstracts from Flybase')
        print('Connecting to Flybase')
        df=pd.io.sql.read_frame(qry,con)
        print('Found %i records'%df.shape[0])
    return df

# pull all gene ids from references with reference id
def pullFBgenerefs(auth=authent_flybase):
    qry='''
    SELECT pub.pub_id, f.uniquename, f.name
    FROM pub, feature_pub fp, feature f, cvterm cv1, cvterm cv2, organism o
    WHERE
    pub.type_id=cv1.cvterm_id AND
    cv1.name='paper' AND
    pub.pub_id=fp.pub_id AND
    fp.feature_id=f.feature_id AND
    f.type_id=cv2.cvterm_id AND
    cv2.name='gene' AND
    f.uniquename LIKE 'FBgn%' AND
    f.is_obsolete='f' AND
    f.organism_id=o.organism_id AND
    o.species='melanogaster'
    '''
    con = psdb.connect(**auth)
    with con:
        print('Getting gene IDs with reference IDs from Flybase')
        print('Connecting to Flybase')
        df=pd.io.sql.read_frame(qry,con)
        print('Found %i records'%df.shape[0])
    return df

'''
Functions for putting flybase data into (local) MySQL database
Note: These recreate entire table from scratch rather than update
'''
# add flybase genes to local database
# note: default MySQL configuration has a size limit to execute lengths
# to work around this, we chunk up into blocks
# alternatively, you can reset this in administrative settings
def putFBgns(df,auth=authent_local,blocsz=10000):
    con = mdb.connect(**auth)
    # add as gene table to local db
    vals= [ tuple([ None if pd.isnull(v) else v for v in rw]) for rw in df.values ]
    sql="""INSERT INTO flygenes(FBid,name)
            VALUES(
            %s,%s)
            """
    nadds=len(vals)*1.0
    nblocs=int(np.ceil(nadds/blocsz))
    print('Adding to flybase genes to local db')
    with con:
        cur=con.cursor()
        
        cur.execute("DROP TABLE IF EXISTS flygenes")
        cur.execute("""CREATE TABLE flygenes(
            geneid INT AUTO_INCREMENT PRIMARY KEY,
            FBid VARCHAR(50), 
            name VARCHAR(100))
            """)       
    #chunk up to keep from bombing
    c=0
    for i in range(nblocs):
        con = mdb.connect(**auth)
        with con:
            cur=con.cursor()
            if (c+blocsz)>nadds:
                cur.executemany(sql,vals[c:])
            else:
                cur.executemany(sql,vals[c:(c+blocsz)])
            c+=blocsz
    print('Done')

# add pub table to local db
def putFBrefs(df,auth=authent_local,blocsz=2000):
    con = mdb.connect(**auth)
    # add as refs table to local db
    vals= [ tuple([ None if pd.isnull(v) else v for v in rw]) for rw in df.values ]
    sql="""INSERT INTO pub(pubid,title,ref,abstract)
            VALUES(
            %s,%s,%s,%s)
            """
    nadds=len(vals)*1.0
    nblocs=int(np.ceil(nadds/blocsz))
    print('Adding to references to local db')
    with con:
        cur=con.cursor()
        
        cur.execute("DROP TABLE IF EXISTS pub")
        cur.execute("""CREATE TABLE pub(
            pubid INT PRIMARY KEY,
            title TEXT, 
            ref TEXT,
            abstract TEXT)
            """)       
    #chunk up to keep from bombing
    c=0
    for i in range(nblocs):
        con = mdb.connect(**auth)
        with con:
            cur=con.cursor()
            if (c+blocsz)>nadds:
                cur.executemany(sql,vals[c:])
            else:
                cur.executemany(sql,vals[c:(c+blocsz)])
            c+=blocsz
    print('Done')

# put gene references into local db
def putFBgenerefs(df,auth=authent_local,blocsz=10000):
    con = mdb.connect(**auth)
    # add as gene references table to local db
    vals= [ tuple([ None if pd.isnull(v) else v for v in rw]) for rw in df.values ]
    sql="""INSERT INTO generefs(pubid,FBid,symb)
            VALUES(
            %s,%s,%s)
            """
    nadds=len(vals)*1.0
    nblocs=int(np.ceil(nadds/blocsz))
    print('Adding gene references to local db')
    with con:
        cur=con.cursor()
        
        cur.execute("DROP TABLE IF EXISTS generefs")
        cur.execute("""CREATE TABLE generefs(
            generef_id INT AUTO_INCREMENT PRIMARY KEY,
            pubid INT,
            FBid VARCHAR(50),
            symb VARCHAR(100))
            """)       
    #chunk up to keep from bombing
    c=0
    for i in range(nblocs):
        con = mdb.connect(**auth)
        with con:
            cur=con.cursor()
            if (c+blocsz)>nadds:
                cur.executemany(sql,vals[c:])
            else:
                cur.executemany(sql,vals[c:(c+blocsz)])
            c+=blocsz
    print('Done')


'''
Main function
Run buildGRNdb to build database
'''
if __name__ == '__main__':
    # build database
    print('Building GeneRefNet database...')
    df=pullFBgns()
    putFBgns(df)
    df=pullFBrefs()
    putFBrefs(df)
    df=pullFBgenerefs()
    putFBgenerefs(df)
    print('Done')
    