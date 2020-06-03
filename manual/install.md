# monolis solid

- monolis structural analysis solver

## インストール

### 0. インストールの準備

Gitlab の [SSH keys](https://gitlab.com/profile/keys) から、公開鍵を登録する。
例えば、[Qiita:【GitLab】SSH認証キー（公開鍵）を登録する](https://qiita.com/CUTBOSS/items/462a2ed28d264aeff7d5) などに詳しい。


### 2. クローン

Gitlab から monolis_solid リポジトリをローカルにクローンする。

```
$ git clone git@gitlab.com:morita/monolis_solid.git
```

monolis_solid ディレクトリに移動する。

```
$ cd monolis_solid
```

### 3. ライブラリのコンパイル

monolis_solid プログラムの線形ソルバは monolis （並列有限要素法向け線形ソルバライブラリ）を利用する。
インストールは、`./install_lib.sh` を実行すれば可能である。
SFEM プログラムは線形ソルバを意識することなく実装されている。
詳細を知りたい場合は、[Gitlab/monolis](https://gitlab.com/morita/monolis) を参照のこと。

```
$ ./install_lib.sh
```

### 4. monolis_solid のコンパイル

`make` コマンドを実行する。

```
$ make
```

インストールが成功していれば、`monolis_solid/bin` に `monolis_solid` が生成されている。
以下のコマンドで確認できる。

```
$ ls bin
```

## サンプル実行

`monolis_solid/example` ディレクトリに以下サンプルが準備されている。

```
circlar_hole # 二次元円孔付き半無限板
```

例えば、以下のようなコマンドで、SFEM プログラムを実行できる。

```
$ cd example/circlar_hole
$ ../../bin/monolis_solid
```
