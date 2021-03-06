---
# Front matter
lang: ru-RU
title: "Отчет по лабораторной работе №4: Модель гармонических коллебаний"
subtitle: "*дисциплина: Математическое моделирование*"
author: "Швец С., НФИбд-03-18"

# Formatting
toc-title: "Содержание"
toc: true # Table of contents
toc_depth: 2
lof: true # List of figures
#lot: true # List of tables
fontsize: 12pt
linestretch: 1.5
papersize: a4paper
documentclass: scrreprt
polyglossia-lang: russian
polyglossia-otherlangs: english
mainfont: PT Serif
romanfont: PT Serif
sansfont: PT Sans
monofont: PT Mono
mainfontoptions: Ligatures=TeX
romanfontoptions: Ligatures=TeX
sansfontoptions: Ligatures=TeX,Scale=MatchLowercase
monofontoptions: Scale=MatchLowercase
indent: true
pdf-engine: luatex
header-includes:
  - \linepenalty=10 # the penalty added to the badness of each line within a paragraph (no associated penalty node) Increasing the value makes tex try to have fewer lines in the paragraph.
  - \interlinepenalty=0 # value of the penalty (node) added after each line of a paragraph.
  - \hyphenpenalty=50 # the penalty for line breaking at an automatically inserted hyphen
  - \exhyphenpenalty=50 # the penalty for line breaking at an explicit hyphen
  - \binoppenalty=700 # the penalty for breaking a line at a binary operator
  - \relpenalty=500 # the penalty for breaking a line at a relation
  - \clubpenalty=150 # extra penalty for breaking after first line of a paragraph
  - \widowpenalty=150 # extra penalty for breaking before last line of a paragraph
  - \displaywidowpenalty=50 # extra penalty for breaking before last line before a display math
  - \brokenpenalty=100 # extra penalty for page breaking after a hyphenated line
  - \predisplaypenalty=10000 # penalty for breaking before a display
  - \postdisplaypenalty=0 # penalty for breaking after a display
  - \floatingpenalty = 20000 # penalty for splitting an insertion (can only be split footnote in standard LaTeX)
  - \raggedbottom # or \flushbottom
  - \usepackage{float} # keep figures where there are in the text
  - \floatplacement{figure}{H} # keep figures where there are in the text
---

# Введение


## Цель
  Изучить и построить модель линейного гармонического осциллятора





# Терминология. Условные обозначения

**Линейный гармони́ческий осцилля́тор**  — система, которая при выведении её из положения равновесия испытывает действие возвращающей силы $F$, пропорциональной смещению $x$


## Уравнение свободных колебаний гармонического осциллятора

Линеное однородное дифференциальное уравнение второго порядка, является примером динамической системы:

$$\ddot{x}=2\gamma\dot{x}+\omega_0^2x$$  
Обозначения:

$\ddot{x} =\frac{\partial^2x}{\partial t^2}$,  $\dot{x} = \frac{\partial x}{\partial t}$

- $x$ -переменная, описывающая состояние системы(смещение грузаб заряд конденсатора и т.д)
- $\gamma$ - характеризует потерю энергии(трение  в  механической системе, сопротивление в контуре)
- $\omega_0$ -собственная частота колебаниц
- $t$ - время

## Уравнение осциллятора при остутствие потерь в системе($\gamma$ =0)

Энергия колебаний такого консервативного  сохраняется во времени:
$$\ddot{x}+\omega_0^2 x =0$$
Зададим начальные условия:
$$\begin{cases}
x(t_0) = x_0
\\
\dot{x}(t_0)=y_0
\end{cases}$$

Уравнение  второго порядка:

$$\begin{cases}
\dot{x}=y
\\
\dot{y} = -\omega_0^2
\end{cases}$$

Начальные условия примут вид:

$$\begin{cases}
x(t_0) = x_0
\\
y(t_0)=y_0
\end{cases}$$


**Фазовая плоскость** - двумерное пространство в которм "движется" решение

**Фазовая траектория** - гладкая кривая в фазовой плоскости - решение уравнения движения как функции времени.

**Фазовый портрет** -  картина, образованная набором фазовых траекторий, когда множество различных решений  можно изобразить на одной фазовой плоскости.



# Выполнение лабораторной работы

## Формулировка задачи:

**Вариант 7**

Постройте фазовый портрет гармонического осциллятора и решение уравнения гармонического осциллятора для следующих случаев:

1. Колебания гармонического осциллятора без затуханий и без действий внешней силы:
$\ddot{x}+7x = 0$

2. Колебания гармонического осциллятора c затуханием и без действий внешней силы:
$\ddot{x}+2\dot{x}+6x = 0$

3. Колебания гармонического осциллятора c затуханием и под действием внешней силы:
$\ddot{x}+5\dot{x}+x = cos(3t)$

На интервале $t \in[0;5]$(шаг 0.05) с начальными условиями $x_0 = -1, y_0 =- 1$

## Решение


Учитывая начальные условия и интервал с шагом $0.05$ построим фазовый портрет гармонического осциллятора и решение уравнения гармонического осциллятора для следующих случаев

1. Колебания гармонического осциллятора без затуханий и без действий внешней силы(рис. -@fig:001)(рис. -@fig:002)

*Решение на  Julia:*  
```

#функция осциляции
function portret(w, g, x0, y0)

    function SDU(du,u,p,t)
        du[1] = u[2]
        du[2] = -w*w*u[1]-g*u[2]-f(t)
    end

    u0 = [x0, y0]
    tspan = (0.0, 25)

    prob = ODEProblem(SDU, u0, tspan)
    sol = solve(prob, RK4(),reltol=1e-6, timeseries_steps = 0.05)

    N = length(sol.u)
    J = length(sol.u[1])

    U = zeros(N, J)

    for i in 1:N, j in 1:J
        U[i,j] = sol.u[i][j]
    end
    U
end
#грфики
f(t) = 0
ans1 = portret(7, 0,1, 1);

Plots.plot(ans1)

set_default_plot_size(30cm, 20cm)
 Gadfly.plot(x = ans1[:,1], y = ans1[:,2],
        Guide.title("Колебания без затухания без действия внешней силы"))
```




![Колебания гармонического осциллятора без затуханий и без действий внешней силы](01.png){ #fig:001 width=70% }

![Решение уравнения для модели гармонического осциллятора без затуханий и без действий внешней силы](011.png){ #fig:002 width=70% }


2. Колебания гармонического осциллятора(рис. -@fig:003) c затуханием и без действий внешней силы(рис. -@fig:004):
$\ddot{x}+13\dot{x}+2x = 0$

Параметры:

  - $\omega = \sqrt{2}$  
  - $\gamma = 13$

*Решение, реаллизованное с помощью  Julia:*  код  аналогичен, меняем лишь вывод
```
set_default_plot_size(30cm, 20cm)
 Gadfly.plot(x = ans1[:,1], y = ans1[:,2],
        Guide.title("Колебания без затухания без действия внешней силы"))

        Plots.plot(ans2)

```


![Колебания гармонического осциллятора c затуханием и без действий внешней силы](02.png){ #fig:003 width=70% }

![Решение уравнения для модели гармонического осциллятора с затуханиями и без действий внешней силы](022.png){ #fig:004 width=70% }


3. Колебания гармонического осциллятора(рис. -@fig:005) c затуханием и под действием внешней силы(рис. -@fig:006):
$\ddot{x}+0.8\dot{x}+1.8x = 2.8sin(8t)$

Параметры:

- $\omega = \sqrt{1.8}$  
- $\gamma = 0.8$

Начальная функция:

- $f(t) = 2.8sin(8t)$



*Решение, реаллизованное с помощью  Julia:*  


```
f(t) = cos(3t)
ans3 = portret(1, 5, -1, 1)
set_default_plot_size(40cm, 20cm)
 Gadfly.plot(x = ans3[:,1], y = ans3[:,2],
        Guide.title("Колебания c затуханием и под действием внешней силы"))

  Plots.plot(ans3)

```


![Колебания гармонического осциллятора c затуханием и действием внешней силы](03.png){ #fig:005 width=70% }

![Решение уравнения для модели гармонического осциллятора с затуханиями и с воздействием внешней силы](033.png){ #fig:006 width=70% }



# Выводы

Мы изучили модель линейного гармонического коллебания и построили ее фазовую траекторю и график решения
