"""
弾性線材(トラス)要素のモデル化について
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"""
###############################################################################
#
# 弾性線材(トラス)要素の剛性方程式は以下になります。
#
# .. math::
#    \begin{bmatrix}\dfrac{EA}{L} & -\dfrac{EA}{L}\\
#    -\dfrac{EA}{L} & \dfrac{EA}{L}
#    \end{bmatrix}\begin{Bmatrix}u_{1}\\
#    u_{2}
#    \end{Bmatrix}=\begin{Bmatrix}F_{1}\\
#    F_{2}
#    \end{Bmatrix}
#
# ここで、 :math:`E` はヤング率、 :math:`A` は断面積、 :math:`L` は要素長さです。
# また、 :math:`u_{1}` と :math:`u_{2}` はそれぞれ節点1と節点2の変位を表し、 :math:`F_{1}` と :math:`F_{2}` はそれぞれ節点力を表します。
#
# .. image:: diagram.png
#
# この定式化が成り立っていることをGetFEMで確認します。

import getfem as gf
import numpy as np

###############################################################################
# パラメータ
# ++++++++++
#
# パラメータの設定は以下の通りに設定します。

E = 205e06  # N/m2
A = 1.0  # m2
L = 10.0  # m

###############################################################################
# メッシュ作成
# ++++++++++++
# :math:`X` 方向に :math:`L` の長さのメッシュを作成します。

mesh = gf.Mesh("cartesian", np.array([0.0, L]))
print(mesh)

###############################################################################
# 有限要素法の定義と積分法
# ++++++++++++++++++++++++
#
# 有限要素法の定義を行います。
# 今回は1次元のため自由度を :math:`1` とします。

elements_degree = 1
mfu = gf.MeshFem(mesh, 1)
mfu.set_classical_fem(elements_degree)

mim = gf.MeshIm(mesh, pow(elements_degree, 1))

###############################################################################
#
# モデル定義
# ++++++++++
#
# 実数モデルを定義し変位(1自由度)を定義します。

md = gf.Model("real")
md.add_fem_variable("u", mfu)

###############################################################################
#
# 線材(トラス要素の定義)
# ++++++++++++++++++++++
#
# 1次元の場合、線材(トラス要素)は次のように定義できます。

md.add_initialized_data("D", [E * A])
md.add_generic_elliptic_brick(mim, "u", "D")


###############################################################################
#
# 剛性行列の確認
# ++++++++++++++
#
# 作成したモデルの接線行列(Tangent Matrix)を確認してみましょう。
# 以下の式の形になっていることを確認してください。
#
# .. math::
#    \begin{bmatrix}\dfrac{EA}{L} & -\dfrac{EA}{L}\\
#    -\dfrac{EA}{L} & \dfrac{EA}{L}
#    \end{bmatrix}
#

md.assembly()
print(md.tangent_matrix().full())

###############################################################################
#
# 斜めの場合のトラスのモデル化について
# ++++++++++++++++++++++++++++++++++++
#
# 次にトラスが斜めの場合のモデル化について確認してみましょう。
# 2次元で斜め60度に回転したモデルを考えます。
# メッシュ作成の際には、まず空の2次元メッシュを作成し、それに要素を点とセルを追加していきます。
#

theta = np.pi / 3.0

mesh = gf.Mesh("empty", 2)
mesh.add_convex(
    gf.GeoTrans("GT_PK(1,1)"), [[0.0, L * np.cos(theta)], [0.0, L * np.sin(theta)]],
)
print(mesh)

###############################################################################
#
# モデル定義
# ++++++++++
#
# 以下の実装は各自由度に1次元の離散化をするのみですので下記のようになります。
#
# .. math::
#    \begin{bmatrix}\dfrac{EA}{L} & 0 & -\dfrac{EA}{L} & 0\\
#    0 & \dfrac{EA}{L} & 0 & -\dfrac{EA}{L}\\
#    -\dfrac{EA}{L} & 0 & \dfrac{EA}{L} & 0\\
#    0 & -\dfrac{EA}{L} & 0 & \dfrac{EA}{L}
#    \end{bmatrix}
#
# ゆえに、回転の影響を考慮することができません。

elements_degree = 1
mfu = gf.MeshFem(mesh, 2)
mfu.set_classical_fem(elements_degree)

mim = gf.MeshIm(mesh, pow(elements_degree, 2))
md = gf.Model("real")
md.add_fem_variable("u", mfu)

md.add_initialized_data("D", [E * A])
md.add_generic_elliptic_brick(mim, "u", "D")

md.assembly("build_matrix")
np.set_printoptions(precision=3)
print(md.tangent_matrix().full())

###############################################################################
#
# 要素座標系から全体剛性回転を次式で定義します。
# 
# .. math::
#    R = \begin{bmatrix}\cos\theta & \sin\theta & 0 & 0\\
#    -\sin\theta & \cos\theta & 0 & 0\\
#    0 & 0 & \cos\theta & \sin\theta\\
#    0 & 0 & -\sin\theta & \cos\theta
#    \end{bmatrix}
#
# 要素座標系の剛性行列を以下で定義します。
#
# .. math::
#    K = \begin{bmatrix}\dfrac{EA}{L} & 0 & -\dfrac{EA}{L} & 0\\
#    0 & 0 & 0 & 0\\
#    -\dfrac{EA}{L} & 0 & \dfrac{EA}{L} & 0\\
#    0 & 0 & 0 & 0
#    \end{bmatrix}
#
# すると剛性行列は以下の式で与えられます。
#
# .. math::
#    R^T K R
#
# 要素座標系の剛性行列は弱定式言語で表現できますが、回転行列をどのように実現するかは考察中です。
