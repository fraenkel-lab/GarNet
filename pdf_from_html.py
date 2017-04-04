from PyQt4.QtGui import QApplication, QPrinter
from PyQt4.QtWebKit import QWebView

def print_pdf(html, destination):
	app = QApplication(sys.argv)

	web = QWebView()
	web.setHtml(html)

	printer = QPrinter()
	printer.setPageSize(QPrinter.A4)
	printer.setOutputFormat(QPrinter.PdfFormat)
	printer.setOutputFileName(destination)
	web.print_(printer)

	app.exit()
