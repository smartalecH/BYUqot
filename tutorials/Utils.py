import tensorflow as tf

linear = lambda x: x

def make_layer(act_func, input_val, input_num, output_num):
    W = tf.Variable(tf.random_normal([input_num, output_num], stddev=0.03), dtype=tf.float32)
    b = tf.Variable(tf.random_normal([output_num], stddev=0.03), dtype=tf.float32)
    layer = act_func(tf.matmul(input_val,W) + b)
    return layer

def make_r(y, pred):
    if y.shape[1] == 2:
        total_error = tf.reduce_sum(tf.square(tf.subtract(y[:,0:1], tf.reduce_mean(y[:,0:1])))) + tf.reduce_sum(tf.square(tf.subtract(y[:,1:2], tf.reduce_mean(y[:,1:2]))))
        unexplained_error = tf.reduce_sum(tf.square(tf.subtract(y[:,0:1], pred[:,0:1]))) + tf.reduce_sum(tf.square(tf.subtract(y[:,1:2], pred[:,1:2])))
        R_squared = tf.subtract(1.0, tf.divide(unexplained_error, total_error))
        return R_squared
    elif y.shape[1] == 1:
        total_error = tf.reduce_sum(tf.square(tf.subtract(y, tf.reduce_mean(y))))
        unexplained_error = tf.reduce_sum(tf.square(tf.subtract(y, pred)))
        R_squared = tf.subtract(1.0, tf.divide(unexplained_error, total_error))
        return R_squared
    else:
        raise TypeError("Weird Shape for R Squared")
