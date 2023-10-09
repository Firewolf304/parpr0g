f [ -f /etc/bashrc ]; then
        . /etc/bashrc
fi
alias treelabs='find . -not -name ".*" | sed -e "s/[^-][^\/]*\// |/g" -e "s/|\([^ ]\)/|-\1/"'
if [ -f ~/.bashrc ]; then
        echo "Tree home page:"
        treelabs
fi
