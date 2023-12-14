import tensorflow as tf


@tf.function
def leq_objective(value, threshold):
    '''
    Realises a "lesser-equal" obective.

    The objective values is zero when the value is below the threshold,
    else it is the distance to the threshold.
    '''
    return tf.maximum(tf.constant(0, dtype='float32'), value - threshold)

@tf.function
def geq_objective(value, threshold):
    '''
    Realises a "greater-equal" obective.

    The objective values is zero when the value is above the threshold,
    else it is the distance to the threshold.
    '''
    return tf.maximum(tf.constant(0, dtype='float32'), threshold - value)

@tf.function
def l1_objective(value, threshold):
    return tf.reduce_sum(tf.abs(value - threshold))

@tf.function
def l2_objective(value, threshold):
    return tf.math.sqrt(tf.reduce_sum((value - threshold)**2))
