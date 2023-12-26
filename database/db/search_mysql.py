import pymysql
config = {
    'host': '127.0.0.1',
    'port': 3306,
    'user': 'root',
    'password': '123456',
    'database': 'database',
    'charset': 'utf8',
}

def search_data(name1, name2, name3, name4):
    db = pymysql.connect(**config)
    cursor = db.cursor()
    sql = "SELECT * FROM " + str(name1) + "_" + str(name2) + "_" + str(name3) + "_" + str(name4)
    cursor.execute(sql)
    data = cursor.fetchall()
    print(data)
    print(type(data))
    db.close()
    return data
