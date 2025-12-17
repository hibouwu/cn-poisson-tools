# cn-poisson-tools
Codes C pour la résolution 1D de l’équation de la chaleur/Poisson (schéma différences finies), avec BLAS/LAPACK, stockage bande, méthodes directes et itératives, plus scripts de validation et mesures de performance.

## Cloner le dépôt / 克隆仓库

```bash
git clone https://github.com/hibouwu/cn-poisson-tools.git
cd cn-poisson-tools
```

## Utilisation avec Docker / 使用 Docker

```bash
# 构建镜像
docker build -f docker/Dockerfile --progress plain -t tp-cn:latest .

# 进入容器编译/运行（将当前目录挂载进去）
sudo docker run --rm -it -v "$(pwd)":/workspace -w /workspace tp-cn:latest make

# 或进入交互 shellDev
sudo docker run -it --rm -v "$(pwd)":/app:z tp-cn:latest /bin/bash
```

Validation de l'installation des bibliothèques BLAS/LAPACK avec le code de test `test_blas_lapack.c` :

```bash
ldconfig -p | grep lapack
```

```txt
liblapacke.so.3 (libc6,x86-64) => /lib/x86_64-linux-gnu/liblapacke.so.3
liblapacke.so (libc6,x86-64) => /lib/x86_64-linux-gnu/liblapacke.so
liblapack.so.3 (libc6,x86-64) => /lib/x86_64-linux-gnu/liblapack.so.3
liblapack.so (libc6,x86-64) => /lib/x86_64-linux-gnu/liblapack.so
```

```bash
# 在仓库根目录进入容器
docker run --rm -it -v "$PWD":/workspace -w /workspace tp-cn:latest bash

# 容器内
mkdir -p bin          # 如果还没建
make                  # 默认目标，生成所有可执行
# 或按需：
# make tp_testenv
# make tpPoisson1D_direct
# make tpPoisson1D_iter

# 运行可执行
./bin/tp_testenv
./bin/tpPoisson1D_direct # 直接法 ./bin/tpPoisson1D_direct [0|1|2]
./bin/tpPoisson1D_iter # 迭代法 ./bin/tpPoisson1D_iter [0|1|2]
```
