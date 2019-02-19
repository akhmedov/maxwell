### UNIT TEST TODO LIST ###

- Math::lommel test
+ numerical derivative test of derivative of I1 and I2
- manual test on maxwell.conf
- caclucaliton manager vs single stream test
+ I1 I3 I4 I5 test on Integral::simpson() method
+ test of getters and setters

### BUG LIST ###

- save mode does not work from plot_model update
+ phi select query in data base
+ mysql connection on Linux; bad_typeid(), race_condition, database NULL fill 
+ config object is not updated if called from main(). Investigate?
- fix noise type configuration by precompiler define-vars
+ computation manager main thread sleep mode (efficiency?)
- numerical linear vs analitical linear
- refactor GnuPlot::write_commants_to_script()
+ refactor of --arg parther in main.cpp
+ refactor and restructur of main.cpp; code is not clean
+ refactor of ready_model.cpp; code copy-paste
+ refactor of manager.cpp; code copy-paste
+ refactor of manager.cpp: online sotring and wait function in get_value()
- extract all cout, cerr and exaptions handling to a separate CLI class
+ catch illegal answers in model
+ Config conflict: Conflict::TEMS_NUM == MissileField::STATIC_TERMS_NUMBER 

- warnings on GNU GCC 5.4.0 20160607
- Makefile refactor; test and maxwell build without rm command
+ implement functionality of Config::print_progress()
+ implement functionality of Config::plot_color_map()
+ implement functionality of Config::plot_grid()
+ implement functionality of Config::plot_baund_cage()
+ implement functionality of Config::plot_device()
+ implement functionality of Config::path_gnuplot_binary()
+ implement functionality of Config::gnp_script_path()
- implement functionality of Config::magnetic_term_num()
+ implement functionality of Config::float_bitrate()

### TODO LIST ###

- stringable enum class
+ secure culculation flag
+ numerics implementaion if Hy(ct,r) in linear medium
- smart pointers all over the code
+ numerical and analitic implementation of E(ct) by CLI option
+ Implement calculation DB
+ Implement logger class
- LIST ALL tested methods next to respount include
+ implement make help option in Makefile
+ gnuplot font size
+ refactor of void Manager::call (AbstractField, string)
+ add THREAD_NUMBER option to maxwell.conf
+ calculation manager progress bar
+ configuration file for global paramiters

+ Optimization of Simpson four-dim algorithm
+ add funvtion MissileEffect::static_field(rho, phi, z, eps = 10e-10)
+ add function GnuPlot::multi_plot2d()
- add function GnuPlot::animation_plot2d()
- add function GnuPlot::animation_plot3d()
+ noise field  AbstuctField::WightAdditiveGaussNoise() 

- call gnuplot from enter point without system() call (version 1.0)
- client server calculation ability through web soket (version 1.0)
- build configuration for Windows 10 (version 1.0)
- implement QT plot device (version 1.0)

### UPDATE LIST FOR c++17 STANDART

- use std::variant in Config:: class for receiver evants
- use std::filsistems::exists() for file check in CLI:: and Manager::

### MYSQL TODO ###

[26.03.2018, 13:34:51] Alex Belinskyi DBA: 1. Какое железо и какой innodb_buffer_pull 
2. Нет индексво на таблицы для WHERE... условий.
[26.03.2018, 13:53:59] rolan a: Привет GL тусовке и спасибо за консультацию))
[26.03.2018, 13:54:06] rolan a: mysql> SELECT @@innodb_buffer_pool_size/1024/1024/1024;
+------------------------------------------+
| @@innodb_buffer_pool_size/1024/1024/1024 |
+------------------------------------------+
|                           0.125000000000 |
+------------------------------------------+
1 row in set (0.00 sec)

если я правильно понял это 124Mb
[26.03.2018, 13:54:43] Alex Belinskyi DBA: 128
[26.03.2018, 13:54:46] Alex Belinskyi DBA: это дефолт
[26.03.2018, 13:54:54] Alex Belinskyi DBA: мало для 100 000
[26.03.2018, 13:55:01] Alex Belinskyi DBA: читает с диска
[26.03.2018, 13:55:06] Alex Belinskyi DBA: индексов нет
[26.03.2018, 13:55:13] Alex Belinskyi DBA: ну кроме как ид
[26.03.2018, 13:55:17] Alex Belinskyi DBA: для джоина
[26.03.2018, 13:55:31] Alex Belinskyi DBA: запросы я в cpp не буду смотреть
[26.03.2018, 13:56:42] Alex Belinskyi DBA: нет времени разбирать что там делается
[26.03.2018, 13:58:29] rolan a: 8Gb хватит для кеша? А про индексы понял, это уже 1/2 ответа. Так что это полезный совет
[26.03.2018, 15:07:32] Alex Belinskyi DBA: сколько памяти на сервере
[26.03.2018, 15:07:38] Alex Belinskyi DBA: мускулю даешь
[26.03.2018, 15:07:42] Alex Belinskyi DBA: 75%
[26.03.2018, 15:07:47] Alex Belinskyi DBA: если там только мускуль
[26.03.2018, 15:14:16] Alex Belinskyi DBA: ну матиматика в помощь. Какой обьем одной строки какая файловая и размер страницы. 
Ну можешь просто,  берешь обьем таблиц их сумарный обьем должен помещатся в кеш. Так как мускуль читает из кеша (innodb_buffer_pooll) в первую очередь.  Большой кеш в 8 Гб то же может тормозить из за количества транзакций. Так что делаются инстансы но опять же если нужно innodb_buffer_pool_instance каждый не более 2 Гб. можно и меньше. Это важно, когда куча паралельных коннектов и т.п.  Но опять же все зависит от того что у тебя в БД происходит и если обьем операционых данных малый.
[26.03.2018, 15:47:41] rolan a: Про индексирование и кеш понял. Результатом доработки похвастаюсь)
[26.03.2018, 15:48:13] Ruslan: результатом доработки - выставься=) Саша вино любит
[26.03.2018, 15:48:56] rolan a: Ваш алкостол был пополнен мной 1.2 раза!))
[26.03.2018, 15:49:35] Alex Belinskyi DBA: ?
[26.03.2018, 15:49:48] Alex Belinskyi DBA: (cwl)
[26.03.2018, 15:49:52] Ruslan: 0,0012 раза?=)
[26.03.2018, 15:50:08] Ruslan: ты не наш пополняй, а Сашин=)
[26.03.2018, 15:50:36] Alex Belinskyi DBA: АТБ -> Baron du vals полусухое) дешево и сердито)
[26.03.2018, 17:49:16] rolan a: 26.03.2018 в 15:49 Ruslan написал (-а):
> 0,0012 раза?=)

порядком-другим ошибся ты)
[26.03.2018, 17:49:31] rolan a: Без проблем, заеду к вам)
[26.03.2018, 21:06:51] rolan a: Саша, скажи какой у тебя стек если не секрет?
[26.03.2018, 21:07:07] Ruslan: Саша
[26.03.2018, 21:07:48] rolan a: Ой, сорри, это все бухлишко))
[26.03.2018, 21:07:56] Ruslan: Что уже пьёшь?
[26.03.2018, 21:08:05] rolan a: Уже нет)
[26.03.2018, 21:08:19] Alex Belinskyi DBA: уже все вылил))
[26.03.2018, 21:09:32] rolan a: на ЗП аспиранта много не купишь)
[26.03.2018, 21:10:18] Ruslan: Сам с синьерской позиции ушёл)
[26.03.2018, 21:10:25] Alex Belinskyi DBA: ты эт не парся выставлятся не надо эт так шутка)
[26.03.2018, 21:10:43] Ruslan: Или не шутка. Неблагодарных людей никто не любит
[26.03.2018, 21:10:48] rolan a: И жалость это плохо)))
[26.03.2018, 21:11:31] Ruslan: Должно произойти нечто ужасное что бы я тебя пожалел(devil)
