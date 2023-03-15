# bomb
对Tinykabomb增加了注释
因为该lab判断交点使用了渐进的方式，每次向光线方向位移一小段距离，根据signed_distance() 和r比较 判断是否相交，所以只要改变signed_distance()就能搞出不同的形状。
```
    Vector3f s = Vector3f(p).normalized() * radius;/*
    float displacement = sin(16 * s.x()) * sin(16 * s.y()) * sin(16 * s.z()) * noise_a*/mplitude;
```
![image](https://user-images.githubusercontent.com/86784362/225347863-01654161-9f37-432c-9925-6386f3aebf09.png)

如果把signed_distance()改为berlin noise 就能实现爆炸的效果。
![image](https://user-images.githubusercontent.com/86784362/225348197-b6fbaa11-1a89-49be-8549-0e3bac6bbd5.png)
