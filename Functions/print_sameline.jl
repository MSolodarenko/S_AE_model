
function println_sameline()
    print("\e[2K") # clear whole line
    print("\e[1G") # move cursor to column 1
    println()
end

function println_sameline(text)
    print("\e[2K") # clear whole line
    print("\e[1G") # move cursor to column 1
    println(text)
end


function print_sameline()
    print("\e[2K") # clear whole line
    print("\e[1G") # move cursor to column 1
    print()
end

function print_sameline(text)
    print("\e[2K") # clear whole line
    print("\e[1G") # move cursor to column 1
    print(text)
end
