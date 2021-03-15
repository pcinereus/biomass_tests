## Usage
##   make -i
DOCS_DIR   := $(addprefix , docs)
#$(info $(DOCS_DIR))
SRC := $(foreach sdir, $(DOCS_DIR), $(wildcard $(sdir)/*.Rmd))
HTML_FILES := $(patsubst %.Rmd, %.html, $(SRC))
#$(info $(SRC))
#$(info $(HTML_FILES))

all: $(HTML_FILES) $(SRC)

$(HTML_FILES): %.html: %.Rmd
	@echo 'Render page'
	@echo $<
	@echo $@
	Rscript -e "library(rmarkdown); render(\"$<\", output_format=\"html_document\", output_options=list(self_contained=TRUE), clean=FALSE)"

clean:
	rm -f *.log *.aux *.md *.out texput.log
