# configure your name and email if you have not done so
git config --global user.email "vladimir.yu.kiselev@gmail.com"
git config --global user.name "wikiselev"

# clone the repository to the book-output directory
git clone -b master \
  https://github.com/wikiselev/bookdown-demo.git \
  book-output
cd book-output
rm -rf *
cp -r ../_book/* ./
git add *
git commit -m "Build the book"
git push origin gh-pages
