CC=gcc -std=c11 -D_GNU_SOURCE
CFLAGS=-c -Wall -Wextra -Werror
TEST_FLAGS=-lcheck -lm -lpthread
GCOV_FLAGS=--coverage -Wall -Wextra -Werror

SOURCES=$(wildcard s21_*.c)
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=s21_matrix

TESTS_SOURCES=test.c
TESTS_OBJECTS=$(TESTS_SOURCES:.c=.o)

ifeq ($(shell uname -s),Linux)
	TEST_FLAGS += -lrt -lsubunit
endif

all: s21_matrix.a

$(EXECUTABLE): $(OBJECTS)
	$(CC) $^ -o $@

.c.o:
	$(CC) $(CFLAGS) $< -o $@

rebuild: clean all

clean:
	rm -rf *.o *.a test tgtests *.gcda *.gcno *.gcov *.html *.css *.out *.info report 

style:
	cp ../materials/linters/.clang-format .clang-format
	clang-format -n *.c *.h
	rm .clang-format

s21_matrix.a: $(OBJECTS)
	ar -rcs s21_matrix.a $(OBJECTS)
	ranlib s21_matrix.a

test: $(TESTS_OBJECTS) s21_matrix.a
	$(CC)  $^ -o test $(TEST_FLAGS)
	./test

valgrind: test
	valgrind --trace-children=yes --track-fds=yes --track-origins=yes --leak-check=full --show-leak-kinds=all -s ./test

gcov_report: $(TESTS_OBJECTS)
	$(CC) $(GCOV_FLAGS) $(TESTS_OBJECTS) $(SOURCES) $(TEST_FLAGS) -o test
	./test
	gcov *.gcno
	lcov -t "test" -o test.info -c -d . --rc lcov_branch_coverage=1
	genhtml -o report test.info --rc lcov_branch_coverage=1
	# gcov s21*.gcda
	gcovr -r . --html --html-details -o report.html
	open report.html
