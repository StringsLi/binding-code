# %%
import tensorflow as tf

# from tensorflow import keras
# from tensorflow.keras import layers
# from tensorflow.python.keras.layers.core import Dropout

# %%
model = tf.keras.models.Sequential([
    tf.keras.Input(shape=(2,)),
    tf.keras.layers.Dense(80, activation='sigmoid'),
    tf.keras.layers.Dense(2)
])

# %%

# x_train=tf.reshape(tf.constant([6.1723,2.3504,7.3475,3.7095,8.5227,3.8935,9.6979,2.9022,10.8731,4.2614,12.0483,4.4453,13.2235,3.454,14.3987,4.8132]),[-1,2])
# y_train=tf.reshape(tf.constant([0,0,1,1,2,2,3,0,4,1,5,2,6,0,7,1]),[-1,2])
# x_test=tf.reshape(tf.constant([15.5739,6.1723]),[-1,2])
# y_test=tf.reshape(tf.constant([8,2]),[-1,2])
import numpy as np

data = np.genfromtxt('NSdata3.csv', delimiter=',')
y = tf.transpose(data[:2])
x = tf.transpose(data[2:])
n = int(len(x) * 0.8)
x_train = x[0:10000]
x_test = x[10000:]
y_train = y[0:10000]
y_test = y[10000:]
# %%
loss_fn = tf.keras.losses.MeanSquaredError()
model.compile(optimizer='adam'
              , loss=loss_fn
              #   ,metrics=['accuracy']
              )
# %%
i = 0
es = 1e-5
maxi = 100
epochs = 20
while True:
    i += 1
    print("#%d..." % i, end="")
    model.fit(x_train, y_train, epochs=epochs, verbose=1)
    print("done")
    loss = model.evaluate(x_train, y_train)
    if loss < es or i >= maxi:
        break

# print("[1,0.1]", model.predict([[0.985, 0.1]]))
# print("[1,0.2]", model.predict([[0.985, 0.2]]))
# print("[1,0.3]", model.predict([[0.985, 0.3]]))
# print("[1,0.4]", model.predict([[0.985, 0.4]]))
# print("[1,0.5]", model.predict([[0.985, 0.5]]))
# print("[1,0.6]", model.predict([[0.985, 0.6]]))
# print("[1,0.7]", model.predict([[0.985, 0.7]]))
# print("[1,0.8]", model.predict([[0.985, 0.8]]))

res = []
for i in range(24, 32):
    x = [1, i * 0.1]
    predict = model.predict([[1, i * 0.1]])
    res.append((x[0], x[1], predict[0][0], predict[0][1]))

import pandas as pd

jj = pd.DataFrame(res)
jj.to_csv('44.csv', index=0)

# %%
# print(model.evaluate(x_test, y_test))
# # %%
# print(model.predict([[1.3259, 3.3203]]))  # 0,2
# # %%
# model.predict([[-1.1374, -1.971],[1.3259, 3.3203],[1.3259, 3.3203]])  # [[-1,-2],[0,2]]
# # %%
# model.predict([[-1.0988, -1.8649]])  # -1.-1.92
# # %%
# model.predict([[0.78459, 2.96976]])  # 0,1.92
#
# # %%
# print(model.predict([[1, 0]]))
# # %%
# model.predict([[-1, 0]])

# %%
# model.save("jy_817.h5")
#
# # %%
# model = tf.keras.models.load_model("jy2.h5")
