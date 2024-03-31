#ifndef TEE_STREAM_CPP
#define TEE_STREAM_CPP

#include <iostream>
#include <fstream>
#include <streambuf>

class TeeStream : public std::streambuf {
public:
    TeeStream(std::streambuf* sb1, std::streambuf* sb2) : sb1(sb1), sb2(sb2) {}

protected:
    // This function is called when outputting each character via the stream.
    virtual int overflow(int c) {
        if (c == EOF) {
            return !EOF;
        } else {
            int const r1 = sb1->sputc(c);
            int const r2 = sb2->sputc(c);
            return (r1 == EOF || r2 == EOF) ? EOF : c;
        }
    }

    // Sync both stream buffers.
    virtual int sync() {
        int const r1 = sb1->pubsync();
        int const r2 = sb2->pubsync();
        return (r1 == 0 && r2 == 0) ? 0 : -1;
    }

private:
    std::streambuf* sb1;
    std::streambuf* sb2;
};

class TeeLogger {
public:
    TeeLogger(const std::string& filepath) : file_stream(filepath),
              tee_buf(std::cout.rdbuf(), file_stream.rdbuf()),
              tee_stream(&tee_buf) {}

    std::ostream& get_stream() {
      return tee_stream;
    }

    std::ofstream& get_file_stream() {
      return file_stream;
    }

private:
    std::ofstream file_stream;
    TeeStream tee_buf;
    std::ostream tee_stream;
};

// Initial code generated using gpt4:
// https://chat.openai.com/c/d51cfec5-e702-4a5b-a471-2487619c4df1
/*
Example usage:

int main() {
    TeeLogger tee_logger("log.txt");
    std::ostream& log_stream = tee_logger.get_stream();

    // This message will go to both std::cout and log.txt
    log_stream << "Hello, world!" << std::endl;
    
    return 0;
}

*/
#endif
