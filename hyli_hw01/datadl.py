import gdown 
tr_path = 'covid.train.csv'  # path to training data
tt_path = 'covid.test.csv'   # path to testing data

url = 'https://drive.google.com/uc?id=19CCyCgJrUxtvgZF53vnctJiOJ23T5mqF'
output = 'covid.train.csv'
gdown.download(url, output, quiet=False)
url = 'https://drive.google.com/uc?id=1CE240jLm2npU-tdz81-oVKEF3T2yfT1O'
output = 'covid.test.csv'
gdown.download(url, output, quiet=False)