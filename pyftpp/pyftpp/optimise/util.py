import numpy as np
# import tensorflow as tf
# import tensorflow_probability as tfp


# @tf.function
# def calculate_dose(dij, w):
#     return tf.sparse.sparse_dense_matmul(dij, w)


def mean_dose(dose, mask):
    '''
    Calculate the mean dose to the given structure.

    The dose and the mask are expected to be 1D tensors of the values
    at each optimisation point and share the same data type.
    '''
    n_voxels = np.sum(mask)
    if n_voxels > 0:
        return np.sum(dose * mask) / n_voxels
    else:
        return 0


# @tf.function
# def maximum_dose(dose, mask):
#     return tf.reduce_max(dose * mask)


# @tf.function
# def interpolate_linearly(x, y, x_eval):
#     '''
#     Parameters
#     ----------
#     x: tf.Tensor
#         x coordinate of the points to be connected by straight lines.
#         Shape: (n,)
#         Assumed to be ordered in ascending order.
#     y: tf.Tensor
#         y coordinate of the points to be connected by straight lines.
#         Shape: (n,)
#     x_eval: tf.Tensor
#         Point at which to evaluate the interpolated curve. Scalar.
#         Assumed to be within the range of the evaluated points,
#         i. e. x_eval < x[-1] and x_eval > x[0].
#     '''
#     if x_eval <= x[0]:
#         return y[0]
#     if x_eval >= x[-1]:
#         return y[-1]
#     # tf.debugging.Assert(
#     #     x_eval > x[0],
#     #     ['evaluation point below interval', x_eval, x[0], x[-1]])
#     # tf.debugging.Assert(
#     #     x_eval < x[-1],
#     #     ['evaluation point above interval', x_eval, x[0], x[-1]])

#     i = 0
#     while x_eval > x[i]:
#         i += 1
    
#     x0 = x[i-1]
#     x1 = x[i]
#     y0 = y[i-1]
#     y1 = y[i]
    
#     epsilon = tf.constant(1e-6, dtype=tf.float32)
#     return y0 + (x_eval - x0) * (y1 - y0) / (x1 - x0 + epsilon)


# @tf.function
# def dvh_d(dose, mask, v):
#     '''
#     Calculates D_<v>%, i. e. the dose to the coldest v% of the voxels in the
#     given structure.
    
#     Parameters
#     ----------
#     v: tf.float32
#         The volume in percent. Must be in [0, 100]. Can be a scalar or
#         of shape (N,).

#     Returns
#     -------
#     tf.float32
#         The dose to the coldest v% of the voxels in the given structure.
#     '''
#     structure_dose = tf.boolean_mask(dose, mask)
#     return tfp.stats.percentile(structure_dose, v)


# def dvh_v(dose, mask, d):
#     '''
#     Calculates an approximation to V_<d>Gy, i. e. the volume that receives
#     at least <d>Gy dose.
#     All float values must have type tf.float32.

#     Parameters
#     ----------
#     dose: tf.Tensor
#         The dose distribution of shape (n,).
#     mask: tf.Tensor
#         The binary mask indicating the structure for which to evaluate V_<d>Gy.
#     d: tf.Tensor
#         The dose in Gy at which to evaluate the DVH.

#     Returns
#     -------
#     tf.float32
#         An approximation to V_<d>Gy [%] by interpolating linearly.
#     '''
#     volume_percentiles = tf.constant(np.linspace(0, 100, 401), dtype=tf.float32)
#     dose_values = dvh_d(dose, mask, volume_percentiles)

#     return interpolate_linearly(dose_values, volume_percentiles[::-1], d)