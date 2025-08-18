function foo(::Broadcast.Broadcasted{BcStyle, Axes, F, Args}) where {BcStyle, Axes, F, Args}
    return BcStyle
end

x = [10, -2, 3]
y = [-3, 1, 19]
bc = Broadcast.Broadcasted(+, (x,y))
res = foo(bc)
print(res)
